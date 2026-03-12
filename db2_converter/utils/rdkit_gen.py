# FROM
# http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
# http://pubs.acs.org/doi/abs/10.1021/ci2004658

import argparse
import multiprocessing
import subprocess
import shutil
import operator
import os
from pathlib import Path
import logging
logger = logging.getLogger("db2_converter")

from rdkit import Chem
from rdkit.Chem import (
    MolFromSmiles,
    MolFromMol2File,
    SanitizeMol,
    SDMolSupplier,
    SDWriter,
    AddHs,
    AssignAtomChiralTagsFromStructure,
    AssignStereochemistryFrom3D,
)
from rdkit.Chem.rdDistGeom import EmbedMolecule, EmbedMultipleConfs, ETKDGv3
from rdkit.Chem.ChemicalForceFields import (
    MMFFGetMoleculeForceField,
    MMFFGetMoleculeProperties,
    UFFGetMoleculeForceField,
)
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

from db2_converter.config import config
from db2_converter.utils.utils import run_external_command
from db2_converter.utils.convert import convert_sdf_to_mol2

RD_NAME = "_Name"
CONF_ENERGY = "Conformer Energy"
FORCEFIELDS = {
    "mmff94": lambda mol, *args, **kwargs: MMFFGetMoleculeForceField(
        mol, MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94"), *args, **kwargs
    ),
    "mmff94s": lambda mol, *args, **kwargs: MMFFGetMoleculeForceField(
        mol, MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s"), *args, **kwargs
    ),
    "uff": UFFGetMoleculeForceField,
}
DEFAULT_FORCEFIELD = "mmff94s"


def convert_by_template(template, mol2block):
    # reorder atoms in new mol2block based on template, to keep different conformations have a same set of atom idxes
    # Besides, to ensure they have the same atom types
    newblock = ""
    newblocklines = list()
    lines = template
    mol2lines = mol2block.split("\n")

    xyz = dict()
    atomflag = False
    for line in mol2lines:  # from mol2block
        if "@<TRIPOS>ATOM" in line:
            atomflag = True
            continue
        if "@<TRIPOS>BOND" in line:
            atomflag = False
            continue
        if atomflag:
            items = line.split()
            atomname = items[0]
            x = "%.4f" % float(items[2])
            y = "%.4f" % float(items[3])
            z = "%.4f" % float(items[4])
            xyz[atomname] = (x, y, z)

    atomflag = False
    for line in lines:  # from template
        if "@<TRIPOS>ATOM" in line:
            atomflag = True
            newblocklines.append(line)
            continue
        if "@<TRIPOS>BOND" in line:
            atomflag = False
            newblocklines.append(line)
            continue
        if atomflag:
            items = line.split()
            atomname = items[0]
            x, y, z = xyz[atomname]
            items[2] = x
            items[3] = y
            items[4] = z
            # newline = '    '.join(items) + '\n'
            newline = f"    {int(items[0]):>3d} {items[1]:<8s} {float(items[2]):>10.4f} {float(items[3]):>10.4f} {float(items[4]):>10.4f} {items[5]:<6s}{items[6]:>6s} {items[7]:<8s}{float(items[8]):>12.4f}\n"
            newblocklines.append(newline)
        else:
            newblocklines.append(line)

    newblock = "".join(newblocklines)
    return newblock


def optimize_single_conformer(
    mol, conf_id, interfragment, max_steps, forcefield=DEFAULT_FORCEFIELD
):
    # use force field in RDKit to minimize and calculate conformational energy
    ignore_frag = not interfragment
    get_forcefield = FORCEFIELDS[forcefield]
    forcefield = get_forcefield(
        mol, confId=conf_id, ignoreInterfragInteractions=ignore_frag
    )
    for _ in range(max_steps):  # Run up to ten steps
        if forcefield.Minimize() == 0:  # If convergence
            break
    energy = forcefield.CalcEnergy()
    return energy


def optimize_conformers(
    mol,
    interfragment=True,
    max_steps=10,
    parallelism=None,
    forcefield=DEFAULT_FORCEFIELD,
):
    # use force field in RDKit to minimize and calculate conformational energy
    conf_energy = {}
    conf_ids = [conf.GetId() for conf in mol.GetConformers()]

    if parallelism is None:
        for conf_id in conf_ids:
            energy = optimize_single_conformer(
                mol,
                conf_id,
                forcefield=forcefield,
                interfragment=interfragment,
                max_steps=max_steps,
            )
            conf_energy[conf_id] = energy
    else:
        args = [
            (mol, conf_id, interfragment, max_steps, forcefield) for conf_id in conf_ids
        ]
        pool = multiprocessing.Pool(processes=parallelism)
        energy = pool.map(optimize_single_conformer, *args)
        conf_energy = dict(list(zip(conf_ids, energy)))

    return conf_energy


def rdk_sample_smi(smiles):
    mol = Chem.MolFromSmiles(smiles)
    params = ETKDGv3()
    params.useSmallRingTorsions = True  # not use small ring torsions
    params.ignoreSmoothingFailures = True  # Prevents crashes in some situations
    mol = Chem.AddHs(mol)
    EmbedMolecule(mol=mol, params=params)
    return mol


def generate_conformers(
    mol,
    add_hydrogens=True,
    rmsd_threshold=2.0,  # Arbitrarily selected
    num_conformers=None,  # None means best guess
    parallelism=None,
    forcefield=DEFAULT_FORCEFIELD,
    skipffcalc=True,
    log=logger,
    max_steps=10,
):  # max_steps for force field minimization
    if add_hydrogens:
        # By default, implicit/explicit hydrogens will all be added
        log.info("Adding implicit hydrogens.")
        mol = AddHs(mol)

    AssignAtomChiralTagsFromStructure(mol)  # can avoid chirality problem?
    AssignStereochemistryFrom3D(mol)

    log.info(
        f"Attempting to generate {num_conformers} conformers at an RMSD cutoff of {rmsd_threshold}."
    )
    params = ETKDGv3()
    params.pruneRmsThresh = rmsd_threshold
    params.useSmallRingTorsions = True  # not use small ring torsions
    params.ignoreSmoothingFailures = True  # Prevents crashes in some situations

    orig_conf_ids = EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
    log.debug(f"Generated {len(orig_conf_ids)} initial conformations.")

    if not max_steps:  # only calculate
        log.info(f"Calculating energies using {forcefield}.")
    else:
        log.info(f"Optimizing and calculating energies using {forcefield}.")
    if skipffcalc:
        log.info(f"Skip forcefield calculation for time saving...")
        return mol, orig_conf_ids
    else:
        conf_energy = optimize_conformers(
            mol,
            interfragment=True,
            parallelism=parallelism,
            forcefield=forcefield,
            max_steps=max_steps,
        )
        sorted_by_energy = sorted(iter(conf_energy.items()), key=operator.itemgetter(1))
        ### pruneRmsThresh above calculates heavy atoms RMSD (After release 2024-03, before it used all atoms)
        ### Below, AlignMol use all-atom RMSD to filter (FF-optimized) conformers.
        ### Notably, symmetry not taken into consideration for RMSD calculation.
        # TODO: Considering adding torsion-based filter like CCDC to reduce time cost.
        log.debug("Filtering similar conformers (considering hydrogens).")
        selected = []
        min_rmsd, max_rmsd = float("inf"), float("-inf")
        for idx, id_energy in enumerate(sorted_by_energy):
            conf_id, energy = id_energy
            keep = True
            for comp_id, other_energy in sorted_by_energy[idx + 1 :]:
                rmsd = AlignMol(mol, mol, prbCid=comp_id, refCid=conf_id)
                if rmsd <= rmsd_threshold:
                    mol.RemoveConformer(conf_id)
                    keep = False
                    break
                else:
                    if rmsd < min_rmsd:
                        min_rmsd = rmsd
                    if rmsd > max_rmsd:
                        max_rmsd = rmsd
            if keep:
                selected.append(id_energy)
        log.debug(
            f"Removed {len(orig_conf_ids)-len(selected)} after post-optimization RMSD filtering"
        )
        log.debug(f"RMSD: min={min_rmsd} max={max_rmsd}")
        return mol, selected


def dump_conformers_mol2(
    mol,
    output,
    params,
    conf_ids=None,
    energies=None,
    renumber=True,
    template=False,
    log=logger,
):
    if conf_ids is None:
        conformers = mol.GetConformers()
    else:
        conformers = (mol.GetConformer(conf_id) for conf_id in conf_ids)

    # Record state of properties that may be overwritten
    original_name = mol.GetProp(RD_NAME)
    if energies is not None and mol.HasProp(CONF_ENERGY):
        original_energy = mol.GetProp(CONF_ENERGY)
    else:
        original_energy = None

    conformer_names = []

    # Render conformers
    tmpdir = "TMP_DIR"
    shutil.rmtree(tmpdir, ignore_errors=True)
    Path(tmpdir).mkdir(exist_ok=False)
    os.chdir(tmpdir)
    if True:
        count = 0
        for idx, conf in enumerate(conformers):
            conf_id = conf.GetId()
            if renumber:
                conf_idx = idx
            else:
                conf_idx = conf_id

            if energies is not None and conf_id in energies:
                energy = energies[conf_id]
                mol.SetProp(CONF_ENERGY, "{0:0.4f}".format(energy))

            conf_name = f"{original_name}_{conf_idx}"
            mol.SetProp(RD_NAME, conf_name)
            conformer_names.append(conf_name)

            try:
                writer = SDWriter(f"{conf_name}.sdf")
                writer.write(mol, confId=idx)
                writer.close()
                convert_sdf_to_mol2(f"{conf_name}.sdf", f"{conf_name}.mol2")
                # convert by template, so no need
                mol2block = open(f"{conf_name}.mol2").read()
                if template:
                    mol2block = convert_by_template(template, mol2block)
                print(mol2block, file=output, end="")
                count += 1
            except:
                continue

        log.debug(f">>> Success rate in rdkit SDF-to-mol2 conversion: {count} / {idx+1}")

    os.chdir("..")
    shutil.rmtree(tmpdir, ignore_errors=False)

    # Reset changes to mol properties
    mol.SetProp(RD_NAME, original_name)
    if original_energy is not None:
        mol.SetProp(CONF_ENERGY, original_energy)
    else:
        mol.ClearProp(CONF_ENERGY)

    return conformer_names


def rdkit_gen(params, log=logger):
    # Load input molecule
    if hasattr(params, "mol2") and params.mol2 is not None:
        mol = MolFromMol2File(params.mol2, sanitize=False)
    elif hasattr(params, "sdf") and params.sdf is not None:
        mol = SDMolSupplier(params.sdf, sanitize=False)[0]
    else:  # SMILES
        smi = Path(params.smi).read_text().split("\n")[0].split()[0]
        mol = MolFromSmiles(smi, sanitize=False)

    try:
        SanitizeMol(mol)
    except Exception as e:
        log.critical(f"Could not sanitize molecule.")
        raise e

    # Assign user-provided name if applicable
    if params.name is not None:
        mol.SetProp(RD_NAME, params.name)
    elif not mol.HasProp(RD_NAME):
        mol.SetProp(RD_NAME, "Ligand")

    # Generate 3D conformers
    embedded, selected = generate_conformers(
        mol,
        add_hydrogens=params.no_hydrogens, # no_hydronges, then add_hydrogens
        rmsd_threshold=params.rmsd_threshold,
        num_conformers=params.num_conformers,
        parallelism=params.parallelism,
        forcefield=params.forcefield,
        skipffcalc=params.skipffcalc,
        max_steps=params.maxstep,
        log=log,
    )
    if not params.skipffcalc:
        log.debug(f"Conformers selected: {len(selected)}")
        log.debug(
            "Energy: min={0:.4f} kcal/mol max={1:.4f} kcal/mol".format(
                selected[0][1], selected[-1][1]
            )
        )

        # Find lowest-energy conformers
        sorted_by_energy = [item[0] for item in selected]
    else:
        sorted_by_energy = selected

    # Render MOL2 file
    template = open(params.template).readlines()
    with open(params.outmol2, "w") as output:
        names = dump_conformers_mol2(
            embedded,
            output,
            params,
            conf_ids=sorted_by_energy,
            renumber=True,
            template=template,
            log=log,
        )
    if not params.skipffcalc:
        for name, (conf_id, energy) in zip(names, selected):
            log.info("\t{0}: {1:0.4f} kcal/mol".format(name, energy))
    os.remove(params.mol2)


# def rdk_enumerate_smi(insmi, tryEmbedding=True):
def rdk_enumerate_smi(insmi, tryEmbedding=False):
    outsmis = []
    m = Chem.MolFromSmiles(insmi)
    opts = StereoEnumerationOptions(tryEmbedding=tryEmbedding)
    isomers = tuple(EnumerateStereoisomers(m, options=opts))  # len(isomers) = 16
    for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
        outsmis.append(smi)
    return outsmis


def smifile_to_sdffile(insmi, outsdf):
    smi = Path(insmi).read_text().split("\n")[0].split()[0]
    name = Path(insmi).read_text().split("\n")[0].split()[1]
    mol = Chem.MolFromSmiles(smi)
    mol = rdk_sample_smi(smi)
    mol.SetProp("_Name", name)
    w = Chem.SDWriter(outsdf)
    w.write(mol)


def main(args, log=logger):
    parser = argparse.ArgumentParser(
        """RDKit-based conformer generation proof-of-concept.
    This program accepts either a mol2 file or a SMILES string and produces a type-fixed mol2 file
    """
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-m", "--mol2", type=str, help="Mol2 file as input to gererate conformers"
    )
    input_group.add_argument("-s", "--smi", type=str, help=".smi file as input")
    input_group.add_argument(
        "--sdf", type=str, help="sdf file as input to generate conformers"
    )

    input_group2 = parser.add_mutually_exclusive_group(required=True)
    input_group2.add_argument("-o", type=str, help="Output mol2 file")

    parser.add_argument("-N", "--name", type=str, default=None, help="Molecule name")
    parser.add_argument(
        "-H",
        "--no-hydrogens",
        action="store_true",
        default=False,
        help="Do NOT explicitly add implicit Hydrogens to conformers [default: %(default)s]",
    )
    parser.add_argument(
        "-r",
        "--rmsd-threshold",
        type=float,
        default=2.0,
        help="Only accept conformers that have an RMSD of at least this value from previously seen conformers [default: %(default)s",
    )
    parser.add_argument(
        "-n",
        "--num-conformers",
        type=int,
        default=None,
        help="Number of conformers to initially generate [default: auto]",
    )
    parser.add_argument(
        "-F",
        "--forcefield",
        type=str,
        default=DEFAULT_FORCEFIELD,
        choices=list(FORCEFIELDS.keys()),
        help="Forcefield to use for optimization [default: %(default)s]",
    )
    parser.add_argument(
        "-P",
        "--parallelism",
        type=int,
        default=None,
        help="Number of processes to use [default: 1]",
    )
    parser.add_argument(
        "-t", "--template", type=str, default="", help="Mol2 file template"
    )
    parser.add_argument(
        "--maxstep", type=int, default="", help="Max minimization steps"
    )
    params = parser.parse_args(args)

    rdkit_gen(params, log=log)
