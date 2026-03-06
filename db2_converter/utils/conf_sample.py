import os
import itertools
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem.Lipinski import NumRotatableBonds, NumAliphaticRings
import logging

logger = logging.getLogger("db2_converter")

from db2_converter.utils.utils import run_external_command
from db2_converter.utils.rdkit_gen import smifile_to_sdffile, rdkit_gen
from db2_converter.config import config


def get_combinations(lst):
    combinations = []
    for r in range(1, len(lst) + 1):
        for subset in itertools.combinations(lst, r):
            combination = "_".join(subset)
            combinations.append(combination)
    return combinations


available_methods = ["conformator", "bcl", "rdkit", "ccdc", "confgenx"]
available_methods = get_combinations(available_methods)


def get_num_confs_for_mol_v1(smi):
    mol = Chem.MolFromSmiles(smi)
    rb = NumAliphaticRings(mol)
    if rb <= 7:
        return 50
    elif 8 <= rb <= 12:
        return 200
    else:
        return 300


def get_num_confs_for_mol_v2(smi):
    mol = Chem.MolFromSmiles(smi)
    rot_bond_num = NumRotatableBonds(mol)
    aliphatic_ring_num = NumAliphaticRings(mol)
    limit_conf_num = (6**rot_bond_num) * (3**aliphatic_ring_num)
    return limit_conf_num


def get_num_confs_for_mol(smi):
    return get_num_confs_for_mol_v2(smi)


class rdk_params:
    def __init__(self, mol2file, name, max_conf=100, rmsd_thres=0.5):
        self.maxstep = 0  # no MMFF opt embedded in rdkit conformer generation
        self.mol2 = "conformer.TMP.fixed.mol2"
        self.template = "conformer.TMP.fixed.mol2"
        self.outmol2 = mol2file
        self.name = name
        self.no_hydrogens = False
        self.num_conformers = max_conf
        self.rmsd_threshold = rmsd_thres
        self.forcefield = "mmff94s"
        self.skipffcalc = True
        self.parallelism = 1


def conf_sample(
    samplopt,
    number,
    mol2file,
    max_conf,
    bcl_option="",
    confgenx_option="",
    RMSDThresh=0.5,
    log=logging,
):
    log = logging.getLogger("conformational sampling")
    UNICON_EXE = config["all"]["UNICON_EXE"]

    if samplopt == "conformator":
        CONF_EXE = config[samplopt]["CONF_EXE"]
        run_external_command(
            f"{CONF_EXE} -i {number}.smi -o {mol2file} -q 2 -n {max_conf} --hydrogens"
        )

    if samplopt in "bcl":
        BCLBASE = config[samplopt]["BCLBASE"]
        BCL = config[samplopt]["BCL"]
        if bcl_option:
            BCL = config[samplopt][bcl_option]
        log.info(f"BCL variant: {BCL}")
        maxIter = 8000  # as stated in BCL paper
        os.environ["PATH"] = f"{BCLBASE}:{os.environ['PATH']}"
        if "LD_LIBRARY_PATH" in os.environ.keys():
            os.environ["LD_LIBRARY_PATH"] = f"{BCLBASE}:{os.environ['LD_LIBRARY_PATH']}"
        else:
            os.environ["LD_LIBRARY_PATH"] = f"{BCLBASE}:"
        smifile_to_sdffile(f"{number}.smi", "conformer.TMP.sdf")
        run_external_command(
            f"""{BCL} -ensemble_filenames conformer.TMP.sdf \
            -generate_3D \
            -top_models {max_conf} \
            -max_iterations {maxIter} \
            -conformation_comparer SymmetryRMSD {RMSDThresh} \
            -conformers_single_file conformer.{number}.sdf \
            -cluster
            """
        )
        run_external_command(f"{UNICON_EXE} -i conformer.{number}.sdf -o {mol2file}", stderr=subprocess.DEVNULL)

    if samplopt == "ccdc":
        CCDC_PYTHON3 = config[samplopt]["CCDC_PYTHON3"]
        ccdc_pyscript = config[samplopt]["CCDC_pyscript"]
        node = config[samplopt]["CCDC_activate_node"]
        smifile_to_sdffile(f"{number}.smi", "conformer.TMP.sdf")
        currpath = os.getcwd()
        ccdc_command = f"ssh {node} \
            {CCDC_PYTHON3} {ccdc_pyscript} \
            --infile {currpath}/conformer.TMP.sdf \
            --outfile {currpath}/conformer.{number}.sdf \
            --max_conf {max_conf} \
            --nthreads 1 \
            --max_unusual_torsions 2"
        # Since we don't have an unlimited license, only one machine in a cluster can be activated simultaneously.
        run_external_command(ccdc_command)
        run_external_command(f"{UNICON_EXE} -i conformer.{number}.sdf -o {mol2file}", stderr=subprocess.DEVNULL)

    if samplopt == "confgenx":
        CONFGENX = config[samplopt]["CONFGENX"]
        SCHUTILS = config[samplopt]["SCHUTILS"]
        if confgenx_option:
            CONFGENX = config[samplopt][confgenx_option]
        log.info(f"ConfGenX variant: {CONFGENX}")
        smifile_to_sdffile(f"{number}.smi", "conformer.TMP.sdf")
        # HH: remove dont_vary_nitrogens option because nitrogen inversion is not taken care of by the enumeration step
        #run_external_command(
        #    f"""{CONFGENX} \
        #    -profile default \
        #    -num_conformers {max_conf} -max_num_conformers {max_conf} \
        #    -auto_increase_threshold 9 -stereo 1 \
        #    -dont_vary_nitrogens \
        #    -j {number} -drop_problematic -no_cleanup -WAIT \
        #    conformer.TMP.sdf
        #    """
        #)
        run_external_command(
            f"""{CONFGENX} \
            -profile default \
            -num_conformers {max_conf} -max_num_conformers {max_conf} \
            -auto_increase_threshold 9 -stereo 1 \
            -j {number} -drop_problematic -no_cleanup -WAIT \
            conformer.TMP.sdf
            """
        )
        # HH: use structconvert to convert maegz to mol2 in one step
        #run_external_command(
        #    f"{SCHUTILS}/structconvert -imae {number}-out.maegz -osd conformer.{number}.sdf"
        #)
        #run_external_command(f"{UNICON_EXE} -i conformer.{number}.sdf -o {mol2file}", stderr=subprocess.DEVNULL)
        run_external_command(
            f"{SCHUTILS}/structconvert -noarom {number}-out.maegz {mol2file}" 
        )

    if samplopt == "rdkit":
        one_rdk_params = rdk_params(
            mol2file=mol2file, name=number, max_conf=max_conf, rmsd_thres=RMSDThresh
        )
        rdkit_gen(one_rdk_params, log=log)


def rdkit_prep(number, mol2file):
    UNICON_EXE = config["all"]["UNICON_EXE"]
    smifile_to_sdffile(f"{number}.smi", "conformer.TMP.sdf")
    run_external_command(f"{UNICON_EXE} -i conformer.TMP.sdf -o conformer.TMP.mol2", stderr=subprocess.DEVNULL)

    shutil.move("conformer.TMP.mol2", "tmp0.mol2")
    shutil.copy("tmp0.mol2", mol2file)
