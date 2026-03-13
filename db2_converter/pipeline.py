import os
import subprocess
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from posebusters import PoseBusters
import logging

logger = logging.getLogger("DB2 generation")


from db2_converter.utils.rdkit_gen import rdk_enumerate_smi
from db2_converter.utils.utils import (
    next_mol2_lines,
    update_mol2block_from_mol,
    check_mol2_smi,
    run_external_command,
    exist_size,
    raise_errlog,
    ATOMTYPE,
)
from db2_converter.utils.convert import (
    convert_sdf_to_mol2,
    restore_dummy_si_in_mol2,
)
from db2_converter.utils.fixmol2 import fixmol2_by_template, fixmol2_and_du
from db2_converter.utils.prepare import prepare, create_namedata
from db2_converter.utils.conf_sample import (
    conf_sample,
    rdkit_prep,
    get_num_confs_for_mol,
)
from db2_converter.amsol.calc_charge_solvation import calc_charge_solvation
from db2_converter.utils.rmsd import RMSDfilter
from db2_converter.utils.match_frags import (
    f_AlignMolConformers,
    embed_blocks_molconformer,
    find_central_rigid,
    mol_to_ring_frags,
)
from db2_converter.mol2db2 import mol2db2
from db2_converter.mol2db2_py3_strain import mol2db2 as mol2db2_38
from db2_converter.strain.Mol2_Strain import mol2_strain
from db2_converter.config import config

# Basic config parameters
UNICON_EXE = config["all"]["UNICON_EXE"]
ANTECHAMBER = config["all"]["ANTECHAMBER"]

# pre-defined parameters
# Don't change if you do not know what they mean.
# EPS = 0.05  # RMSD likelihood threshold
EPS = 0.15
EPS_modify = 1e-4  # lower kept, upper modified
max_amsol_attempts = 10  # amsol attemp limit
RMSthres = 0.5  # RMSD threshold
sb_smarts = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"
minfragsize = 3


def max_conf_assign(smi, max_conf, max_conf_noring):
    mol = Chem.MolFromSmiles(smi)
    if not len(mol_to_ring_frags(mol, cutsmarts=sb_smarts)):
        max_conf = max_conf_noring
    return max_conf


def sample_tp_unicon(infile,outfile=""):
    if not outfile:
        outfile = infile.stem + ".uni.smi"
        # only generate one smiles
        run_external_command(f"{UNICON_EXE} -i {infile} -o {outfile} -t single -p single")
    return Path(outfile)


def write_enumerated_smifile(inlines, outfile, method):
    logger.info(
        "####### Enumerating possible undefined stereochemistry of input SMILES... #######"
    )
    newlines = []
    toomany = False
    for line in inlines:
        try:
            smi, name = line.split()[0], line.split()[1]
            outsmis = rdk_enumerate_smi(smi)
            if len(outsmis) == 1:
                newlines.append(line)
            else:
                if len(outsmis) <= 2**5:
                    for i, outsmi in enumerate(outsmis):
                        newlines.append(f"{outsmi} {name}.{i}")
                        logger.info(f">>> {name}.{i} {outsmi}")
                else:
                    logger.error(
                        f">>> SMILES of {name} has too many stereoisomers, so skipped!"
                    )
                    toomany = True
                    faillist = [smi, name, method, "0enumerate_many"]
        except Exception as e:
            logger.error(e)
            logger.error(f">>> SMILES of {name} cannot be parsed")
    logger.info("####### Enumeration finished. #######\n")

    with open(outfile, "w") as f:
        for line in newlines:
            f.write(f"{line}\n")
    if toomany:
        return faillist


def fixmol2_wrapper(
    inmol2file,
    outmol2file,
    templatemol2file="",
    smiles="",
    samplopt="",
    covalent=False,
):
    tmp0mol2 = "tmp0.mol2"
    tmp0mol2_ante = "tmp0.ante.mol2"
    tmp0fixmol2 = "tmp0.mol2.fixed.mol2"
    TMPmol2 = "conformer.TMP.fixed.mol2"
    for i, mol2content in enumerate(next_mol2_lines(inmol2file)):
        with open(f"tmp{i}.mol2", "w") as f:
            f.write("".join(mol2content))
    si_atoms = {}
    antechamber_input = tmp0mol2
    if covalent:
        si_atoms = _swap_dummy_si_for_antechamber(tmp0mol2, tmp0mol2_ante)
        if si_atoms:
            antechamber_input = tmp0mol2_ante
    antechamber_cmd = (
        f"{ANTECHAMBER} -i {antechamber_input} -fi mol2 "
        f"-o {tmp0fixmol2} -fo mol2 -at sybyl -pf y"
    )
    logger.info(">>> Running antechamber: %s", antechamber_cmd)
    result = run_external_command(antechamber_cmd)
    restored_tmp0fixmol2 = tmp0fixmol2
    if exist_size(tmp0fixmol2):
        antechamber_out = f"{tmp0fixmol2}.ante"
        shutil.copy(tmp0fixmol2, antechamber_out)
        if covalent and si_atoms:
            restored_tmp0fixmol2 = f"{tmp0fixmol2}.restored"
            shutil.copy(tmp0fixmol2, restored_tmp0fixmol2)
            restore_dummy_si_in_mol2(restored_tmp0fixmol2, si_atoms)
    logger.info(
        ">>> Antechamber exit code: %s; output exists: %s",
        result.returncode,
        exist_size(tmp0fixmol2),
    )
    if not templatemol2file:
        if not smiles:
            assert smiles != "", "SMILES is not given to fixmol2"
        if exist_size(restored_tmp0fixmol2):
            # Apply fixes before SMILES validation so nitro/charge cleanup can succeed.
            fixmol2_and_du(smiles, restored_tmp0fixmol2)
        check_ok = True
        if not exist_size(restored_tmp0fixmol2):
            check_ok = False
        else:
            check_ok = check_mol2_smi(restored_tmp0fixmol2, smiles)
        logger.info(">>> Antechamber mol2 SMILES check: %s", check_ok)
        if not check_ok:
            # antechamber cannot deal with "[N-]" due to unexpected valence.
            logger.info(">>> Falling back to raw tmp0.mol2 as template")
            shutil.copy(tmp0mol2, tmp0fixmol2)
            restored_tmp0fixmol2 = tmp0fixmol2
            fixmol2_and_du(smiles, restored_tmp0fixmol2)
        templatemol2file = restored_tmp0fixmol2
    if samplopt == "rdkit":
        shutil.move(tmp0fixmol2, TMPmol2)
    else:
        for tmpmol2 in sorted(Path(".").glob("tmp*.mol2")):
            fixmol2_by_template(tmpmol2, templatemol2file)
        # Keep template file for debugging.
        # subprocess.run(f"rm {tmp0fixmol2}", shell=True)
        subprocess.run(f"cat tmp*.mol2 > {outmol2file}", shell=True)
        # Keep tmp* files for debugging.
        # subprocess.run("rm tmp*", shell=True)
    return TMPmol2 if samplopt == "rdkit" else templatemol2file


def _swap_dummy_si_for_antechamber(inmol2, outmol2):
    si_atoms = {}
    mol2blocks = [x for x in next_mol2_lines(inmol2)]
    if not mol2blocks:
        shutil.copy(inmol2, outmol2)
        return si_atoms
    new_blocks = []
    for mol2block in mol2blocks:
        atom_start = mol2block.index("@<TRIPOS>ATOM\n")
        bond_start = mol2block.index("@<TRIPOS>BOND\n")
        startpart = mol2block[: atom_start + 1]
        atompart = mol2block[atom_start + 1 : bond_start]
        bondpart = mol2block[bond_start:]
        new_atompart = []
        for line in atompart:
            items = line.strip().split()
            if len(items) < 9:
                new_atompart.append(line)
                continue
            atom_id = int(items[0])
            atom_name = items[1]
            atom_type = items[5]
            if atom_type == "Si" or atom_type.startswith("Si"):
                si_atoms[atom_id] = atom_name
                atom_name = f"C{atom_id}"
                atom_type = "C.3"
                new_atompart.append(
                    ATOMTYPE.format(
                        atom_id,
                        atom_name,
                        float(items[2]),
                        float(items[3]),
                        float(items[4]),
                        atom_type,
                        items[6],
                        items[7],
                        float(items[8]),
                    )
                )
            else:
                new_atompart.append(line)
        new_blocks.append("".join(startpart + new_atompart + bondpart))
    with open(outmol2, "w") as f:
        f.write("".join(new_blocks))
    return si_atoms


def chemistrycheck(insmi, inmol2, outmol2, checkstereo=True):
    canonical_smiles = Chem.MolToSmiles(
        Chem.MolFromSmiles(insmi), isomericSmiles=checkstereo
    )
    all_blocks = [x for x in next_mol2_lines(inmol2)]
    all_mols = [Chem.MolFromMol2Block("".join(x), removeHs=True) for x in all_blocks]
    correct_indexes = []
    for i, mol in enumerate(all_mols):
        if checkstereo:
            # Chem.AssignAtomChiralTagsFromStructure(mol)
            """From rdkit:
            NOTE that this function does not check if atoms are chiral centers (i.e. all substituents are different),
            it merely sets the chiral type flags based on the coordinates and atom ordering.
            Use AssignStereochemistryFrom3D() if you want chiral flags only on actual stereocenters.
            """
            Chem.AssignStereochemistryFrom3D(mol)
        smi = Chem.MolToSmiles(mol, isomericSmiles=checkstereo)
        if smi == canonical_smiles:
            correct_indexes.append(i)
        else:
            logger.warning(f"Not consistent:\ngen {smi}\nref {canonical_smiles}")
    logger.info(f">>> {len(correct_indexes)} / {len(all_mols)} passed chemistrycheck.")
    with open(outmol2, "w") as f:
        for i in correct_indexes:
            f.write("".join(all_blocks[i]))


def PB_filter(inmol2, outmol2):
    all_blocks = [x for x in next_mol2_lines(inmol2)]
    run_external_command(
        f"{UNICON_EXE} -i {inmol2} -o {inmol2}.unipb.sdf",
    )
    buster = PoseBusters(config="mol")
    df = buster.bust([f"{inmol2}.unipb.sdf"], full_report=True)
    # df = buster.bust([f"{inmol2}"], full_report=True)
    df["PBvalid_conf"] = (
        df["all_atoms_connected"]
        & df["sanitization"]
        & df["bond_lengths"]
        & df["bond_angles"]
        & df["internal_steric_clash"]
        & df["aromatic_ring_flatness"]
        & df["double_bond_flatness"]
        & df["internal_energy"]
    )
    df = df.reset_index(drop=False)
    correct_indexes = df.index[df["PBvalid_conf"]].tolist()
    logger.info(
        f">>> {len(correct_indexes)} / {len(all_blocks)} passed PoseBusters filter."
    )
    with open(outmol2, "w") as f:
        for i in correct_indexes:
            f.write("".join(all_blocks[i]))


def mmffopt(inmol2, outmol2):
    mol2blocks = [i for i in next_mol2_lines(inmol2) if i]
    newmol2content = []
    for mol2block in mol2blocks:
        newmol = Chem.MolFromMol2Block("".join(mol2block), removeHs=False)
        MMFFOptimizeMolecule(newmol)
        newmol2block = update_mol2block_from_mol(mol2block, newmol)
        newmol2content.append("".join(newmol2block))
    with open(outmol2, "w") as f:
        f.write("".join(newmol2content))


def prepare_mol2(name, smiles, inmol2, mergeiso=True):
    """
    # prapare.py
    ## generate a file called name.txt, one example content is:
    ## name.txt 1 CHEMBL85549 OC[C@H]1OC(N2C=NC3=C(NC4CCCCC4)N=CN=C32)[C@H](O)[C@@H]1O | NO_LONG_NAME
    ## This file will be used as default -n option in mol2db2.py
    """
    if mergeiso:
        name = name.split(".")[0]
    namedata = create_namedata(name=name, smiles=smiles, longname="NO_LONG_NAME")
    # for DUD-E / DUDE-Z when multiple protomers are considered as a single molecule
    # namedata = create_namedata(name=name.split(".")[0].split("_")[0], smiles=smiles, longname="NO_LONG_NAME")
    prepare(inmol2, namedata=namedata, namepath=None, tempdir=".")


def match_and_convert_mol2(
    mol2file,
    extra_fragsindex=[],
    extra_fragsmarts="",
    onlyextrafrags=False,
    reseth=True,
    rotateh=True,
    prefix="output",
    dock38=False,
    covalent=False,
    templatemol2file="",
):
    all_blocks = [x for x in next_mol2_lines(mol2file)]
    mol = Chem.MolFromMol2Block("".join(all_blocks[0]), removeHs=False)
    # Covalent mode uses a Si-centered rigid fragment.
    if covalent:
        si_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 14]
        if si_atoms:
            si_idx = si_atoms[0].GetIdx()
            visited = {si_idx}
            frontier = {si_idx}
            for _ in range(2):
                next_frontier = set()
                for idx in frontier:
                    atom = mol.GetAtomWithIdx(idx)
                    for neighbor in atom.GetNeighbors():
                        nidx = neighbor.GetIdx()
                        if nidx not in visited:
                            visited.add(nidx)
                            next_frontier.add(nidx)
                frontier = next_frontier
            fragsindex = [
                [idx for idx in sorted(visited) if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1]
            ]
        else:
            logger.warning(">>> Covalent mode enabled but no Si atom found; falling back to ring-based fragments.")
            fragsindex = []
    else:
        fragsindex = []
    # Get atom name to serial number mapping
    index_of_name = dict()
    for atom in mol.GetAtoms():
        name = atom.GetPropsAsDict()["_TriposAtomName"]  # extract O1, C2
        index_of_name[name] = atom.GetIdx()
    # Get break down fragments
    if not fragsindex:
        suitable_ring_frags = mol_to_ring_frags(
            mol=mol, cutsmarts=sb_smarts, minfragsize=3
        )
        for frag in suitable_ring_frags:
            findex = []
            for atom in frag.GetAtoms():
                if atom.GetAtomicNum() > 1:  # Heavy atoms
                    name = atom.GetPropsAsDict()["_TriposAtomName"]
                    findex.append(index_of_name[name])
            # Double check to ensure at least 3 heavy atoms in one fragment
            if len(findex) >= minfragsize:
                fragsindex.append(findex)
    # extra rigid fragment #
    if not covalent:
        if extra_fragsindex:
            extra_fragsindex = [extra_fragsindex]
        if extra_fragsmarts:
            extra_fragsindex = []
            for extra_fragindex in mol.GetSubstructMatches(
                Chem.MolFromSmiles(extra_fragsmarts)
            ):
                extra_fragsindex += [list(extra_fragindex)]
        if extra_fragsindex:  # extra_fragsmarts has higher priority than extra_fragsindex
            fragsindex += extra_fragsindex
            if onlyextrafrags:
                fragsindex = extra_fragsindex
    ########################
    if not fragsindex:
        fragsindex = [find_central_rigid(mol)]
        logger.info(f">>> No rigid part can be kept after bond cleavage.")
        logger.info(
            f">>>>>> Will use central atom and its 1st neighbors {fragsindex}..."
        )
    logger.debug(f"fragsindex: {fragsindex}")

    mol_withconfs = embed_blocks_molconformer(all_blocks, removeHs=False)
    # Aligning every frag
    i = 0
    for dirname in ["sdf", "mol2", "db2"]:
        Path(dirname).mkdir(exist_ok=True)
    for index in fragsindex:
        all_index = f_AlignMolConformers(
            mol_withconfs,
            atomIds=index,
            maxIters=500,
            reflect=False,
            eps=EPS,
            eps_modify=EPS_modify,
        )  # maxIters and reflect are for subfunc AlignMolConformers
        # all_index is a list of lists, each list with len>1 means a cluster
        logger.debug("Rigid body index: %s", index)
        logger.debug("Conformers used id: %s", all_index)
        for conf_index in all_index:
            writer = Chem.SDWriter(f"sdf/{prefix}.{i}.sdf")
            for conf in conf_index:
                writer.write(
                    mol_withconfs, confId=conf
                )  # write different aligned confs into separated sdf files
            writer.close()
            convert_sdf_to_mol2(
                f"sdf/{prefix}.{i}.sdf",
                f"mol2/{prefix}.{i}.mol2",
            )
            if templatemol2file and exist_size(templatemol2file):
                fixmol2_by_template(f"mol2/{prefix}.{i}.mol2", templatemol2file)
            # # if output mol2 has format issue, put fixmol2_wrapper here
            # fixmol2_wrapper(
            #     f"mol2/{prefix}.{i}.mol2",
            #     f"mol2/{prefix}.{i}.fixed.mol2",
            #     templatemol2file=mol2file
            # )
            if not dock38:
                ############## DOCK37 db2 from mol2 ##############
                mol2db2.mol2db2_main(
                    mol2file=f"mol2/{prefix}.{i}.mol2",
                    solvfile=f"{prefix}.solv",
                    namefile="name.txt",
                    db2gzfile=f"db2/{prefix}.{i}.db2.gz",
                    timeit=True,
                    reseth=reseth,
                    rotateh=rotateh,
                )
                ##################################################
            else:
                ############## DOCK37 db2 from mol2 ##############
                fixmol2_wrapper(
                    f"mol2/{prefix}.{i}.mol2",
                    f"mol2/{prefix}.{i}.fixed.mol2",
                    templatemol2file=mol2file
                )
                db2in_standard = mol2_strain(f"mol2/{prefix}.{i}.fixed.mol2")
                mol2db2_38.mol2db2_main(
                    db2in_standard,
                    solvfile=f"{prefix}.solv",
                    db2gzfile=f"db2/{prefix}.{i}.db2.gz",
                    timeit=True,
                    reseth=reseth,
                    rotateh=rotateh,
                )
            i += 1
    return i


def gen_conf(
    zinc,
    max_conf,
    max_conf_noring,
    samplopt,
    checkstereo,
    PBfilter=False,
    MMFFopt=False,
    cluster=False,
    limitconf=False,
    extra_fragsindex=[],
    extra_fragsmarts="",
    onlyextrafrags=False,
    keep_max_conf=False,
    reseth=True,
    rotateh=True,
    prefix="output",
    faillist=[],
    bcl_option="",
    confgenx_option="",
    mergeiso=True,
    dock38=False,
    **kwargs
):
    logger.info(f"############### Now dealing with {zinc}... ###############")
    Path(zinc).mkdir(exist_ok=True)
    shutil.copy(f"{zinc}.smi", zinc)
    os.chdir(zinc)
    smi = Path(f"{zinc}.smi").read_text().split("\n")[0].split()[0]
    ringmol = Chem.MolFromSmiles(smi)
    N_ring = len(mol_to_ring_frags(ringmol, cutsmarts=sb_smarts))
    if rotateh:
        NrotHs, multiplier = mol2db2.mol2db2_to_numhyds(f"{zinc}.smi")
    if not keep_max_conf:
        if not N_ring:
            max_conf = max_conf_noring
        if N_ring and rotateh:  # if rotateh, we need to truncate the ensemble size
            if NrotHs in [4, 5]:
                max_conf = max_conf // 30
            if NrotHs in [2, 3]:
                max_conf = max_conf // 3
        if limitconf:
            max_conf = min(get_num_confs_for_mol(smi), max_conf)
    # if keep_max_conf, the max_conf you defined is the actual max_conf for generation
    N_max_conf_in = max_conf

    samplname = samplopt
    samplopts = samplopt.split("_")
    if len(samplopts) > 1:
        logger.info(f"#### {samplname} ensemble sampling... ####")

    templatemol2file = ""
    for samplopt in samplopts:
        mol2file = f"conformer.{zinc}.{samplopt}.mol2"
        fixed_mol2file = f"conformer.{zinc}.{samplopt}.fixed.mol2"

        logger.info(
            f">>> {samplopt} used to sample {max_conf} conformers, checkstereo {checkstereo}, PBfilter {PBfilter}, MMFFopt {MMFFopt}, cluster {cluster}"
        )
        if samplopt == "rdkit":
            rdkit_prep(zinc, mol2file)
        else:
            conf_sample(
                samplopt,
                zinc,
                mol2file,
                max_conf,
                bcl_option,
                confgenx_option,
                log=logger,
            )

        if not exist_size(mol2file):
            error = "1generate"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            continue

        try:
            templatemol2file = fixmol2_wrapper(
                mol2file,
                fixed_mol2file,
                "",
                smi,
                samplopt,
                covalent=kwargs.get("covalent", False),
            )
        except Exception as e:
            logger.error(e)
            error = "2fixmol2"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            continue

        if samplopt == "rdkit":
            try:
                conf_sample(samplopt, zinc, fixed_mol2file, max_conf, log=logger)
            except Exception as e:
                logger.error(e)
                error = "1generate"
                faillist.append([smi, zinc, samplopt, error])
                raise_errlog(error, logger)
                continue

        if not exist_size(fixed_mol2file):
            error = "2fixmol2"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            continue

        ## Chemistry check and filter (Through RDKit)
        if checkstereo:
            logger.info(">>> Chemistry checking with stereochemistry...")
        else:
            logger.info(">>> Chemistry checking without stereochemistry... Be careful!!!")
        try:
            chemistrycheck(
                insmi=smi,
                inmol2=fixed_mol2file,
                outmol2=f"conformer.{zinc}.fixed.filter.mol2",
                checkstereo=checkstereo,
            )
            shutil.move(f"conformer.{zinc}.fixed.filter.mol2", f"{fixed_mol2file}")
        except Exception as e:
            with open(fixed_mol2file, "w") as f:
                pass  # check failed, so we do not want to keep
        if os.path.getsize(fixed_mol2file) == 0:
            error = "3chemistrycheck"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            continue
        else:
            logger.info(">>> Chemistry Checked.")

        ## PoseBusters filter
        if PBfilter:
            logger.info(">>> PoseBusters filter implausible conformers...")
            PB_filter(inmol2=fixed_mol2file, outmol2=f"conformer.{zinc}.fixed.pb.mol2")
            shutil.move(f"conformer.{zinc}.fixed.pb.mol2", f"{fixed_mol2file}")
            logger.info(">>> PoseBusters filter Finished.")
        if os.path.getsize(fixed_mol2file) == 0:
            error = "4PBfilter"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            continue

        ## MMFF optimization
        if MMFFopt:
            logger.info(">>> MMFF optimization...")
            mmffopt(
                inmol2=fixed_mol2file, outmol2=f"conformer.{zinc}.fixed.mmffopt.mol2"
            )
            shutil.move(f"conformer.{zinc}.fixed.mmffopt.mol2", f"{fixed_mol2file}")
            logger.info(">>> MMFF optimization Finished.")

        ## Iterate to get one normal amsol output
        logger.info(">>> AMSOL calculation for atom partial charges and desolvation...")
        if not exist_size(f"{prefix}.solv"):  # keep for ensemble
            mol2lines_list = list(next_mol2_lines(fixed_mol2file))
            for i, mol2lines in enumerate(mol2lines_list):
                logger.debug(f"amsol trial {i}")
                with open(f"{prefix}{zinc}.mol2", "w") as f:
                    f.write("".join(mol2lines))

                prepare_mol2(
                    name=zinc,
                    smiles=smi,
                    inmol2=f"{prefix}{zinc}.mol2",
                    mergeiso=mergeiso,
                )

                # AMSOL ( very complicated in it, but not much work )
                ## with a lot of output log info
                ## ::input: output$number.mol2 (single mol2), total charge used by AMSOL is directly summed read from the APC in it
                ## ::output: {prefix}.solv (amsol output file)
                ## ::output: {prefix}.mol2 (single mol2, w/ charge)
                ## ::output: {prefix}$number.mol2 (single mol2, format slightly changed, w/o charge)
                try:
                    calc_charge_solvation(f"{prefix}{zinc}.mol2")
                except:
                    pass
                if exist_size(f"{prefix}.solv") or i + 1 == max_amsol_attempts:
                    logger.info(">>> AMSOL calculation Finished.")
                    break

            if not exist_size(f"{prefix}.solv"):
                error = "5amsolfail"
                faillist.append([smi, zinc, samplopt, error])
                raise_errlog(error, logger)
                continue
        faillist.append([])
        # If success, faillist will be [], this could be helpful for combinatorial generation

    # collect mol2, rmsd filter
    fixed_mol2file = f"conformer.{zinc}.fixed.mol2"
    fail_name_count = len([onefail for onefail in faillist if onefail])
    if fail_name_count == len(samplopts):
        return faillist, N_max_conf_in, "", "", N_ring, ""
    else:
        origin_mol2_blocks = []
        for samplopt in samplopts:
            if exist_size(f"conformer.{zinc}.{samplopt}.fixed.mol2"):
                mol2_blocks = list(
                    next_mol2_lines(f"conformer.{zinc}.{samplopt}.fixed.mol2")
                )
                origin_mol2_blocks += mol2_blocks
        if not cluster:
            out_mol2_blocks = origin_mol2_blocks
        else:
            logger.info(f">>> RMSD-based Clustering at {RMSthres} Angstrom...")
            out_mol2_blocks = RMSDfilter(zinc, samplopts, RMSthres)
            logger.info(f">>> {len(out_mol2_blocks)} / {len(origin_mol2_blocks)} kept.")
            logger.info(">>> RMSD Clustering finished...")
        with open(fixed_mol2file, "w") as f:
            f.write("".join(["".join(mol2_block) for mol2_block in out_mol2_blocks]))
    N_act_conf = len(out_mol2_blocks)
    if rotateh:
        N_act_conf_out_rotH = N_act_conf * multiplier
    else:
        N_act_conf_out_rotH = N_act_conf

    # Match and convert
    ## ::input: conformer.$number.fixed.mol2 (infile)(multiple conformations)
    ## ::output: dir sdf with separated cluster sdf file
    ## ::output: dir mol with separated cluster mol2 file
    ## ::output: dir db2 with separated db2.gz file
    shutil.rmtree("db2", ignore_errors=True)
    template_for_align = f"{prefix}.mol2"
    db2part_count = match_and_convert_mol2(
        mol2file=fixed_mol2file,
        extra_fragsindex=extra_fragsindex,
        extra_fragsmarts=extra_fragsmarts,
        onlyextrafrags=onlyextrafrags,
        reseth=reseth,
        rotateh=rotateh,
        prefix="output",
        dock38=dock38,
        covalent=kwargs.get("covalent", False),
        templatemol2file=template_for_align,
    )
    # collect
    if Path("db2").exists():
        subprocess.run("cat db2/*.db2.gz > all.db2.gz", shell=True)
    if os.path.getsize("all.db2.gz") != 28:
        logger.info(f">>> {db2part_count} db2 blocks were generated and saved...")
        logger.info(f"############### Finished with {zinc}... ###############\n")
    else:
        error = "9nulldb2gz"
        faillist.append([smi, zinc, samplopt, error])
        raise_errlog(error, logger)

    return (
        faillist,
        N_max_conf_in,
        N_act_conf,
        N_act_conf_out_rotH,
        N_ring,
        db2part_count,
    )
