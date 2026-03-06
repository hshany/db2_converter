import os
import sys
import logging
import subprocess
from pathlib import Path
from db2_converter.config import config
from db2_converter.utils.utils import exist_size, run_external_command
from db2_converter.amsol.make_amsol71_input import make_amsol71_input
from db2_converter.amsol.process_amsol_mol2 import process_amsol_mol2

logger = logging.getLogger("amsol")


# HH: reduce amsoltimelimit from 60 to 10
def calc_charge_solvation(
    mol2file, amsoltimelimit=10, logger=logger
):  # default 1min for each amsol calculation
    AMSOLEXE = config["all"]["AMSOLEXE"]
    OBABELEXE = config["all"]["BABEL_EXE"]

    # Obtain the name of the protonated molecule from the mol2-file name
    constituents = str(mol2file).split(".")
    protonated_molecule_name = f"{constituents[0]}.{constituents[1]}"
    # logger.debug(f"ProtonatedMoleculeName: {protonated_molecule_name}")
    molecule_name = protonated_molecule_name

    logger.debug(
        f"Preparing AMSOL7.1 input for {mol2file} (first: transformation to ZmatMOPAC by openbabel)"
    )
    logger.debug("are correctly set in your environment variables.")

    if exist_size("temp.mol2"):
        logger.warn("Warning: temp.mol2 exists. Rewriting link.")
        # Link the mol2 file to temp.mol2
        os.remove("temp.mol2")
    os.symlink(mol2file, "temp.mol2")

    # Obtain the entire current path
    current_path = os.getcwd()
    logger.debug(f"current_path : {current_path}")

    temp_file = f"{current_path}/temp.mol2"

    # Convert mol2 to ZmatMOPAC format using obabel
    logger.debug(
        f"Running: {OBABELEXE} -i mol2 {temp_file} -o mopin -O {current_path}/temp.ZmatMOPAC"
    )
    run_external_command(
        f"{OBABELEXE} -i mol2 {temp_file} -o mopin -O {current_path}/temp.ZmatMOPAC",
        stderr=subprocess.STDOUT
    )

    # Set PYTHONPATH environment variable
    # Create AMSOL7.1 input files
    # TODO: python as submodule.
    logger.debug(
        f"Make amsol71_input from {current_path}/temp.ZmatMOPAC {molecule_name}"
    )
    make_amsol71_input(f"{current_path}/temp.ZmatMOPAC", molecule_name)

    # Run AMSOL7.1 calculations
    logger.debug("Running AMSOL7.1: SM5.42R (in water solvent)")
    run_external_command(
        f"{AMSOLEXE} < temp.in-wat > temp.o-wat", timeout=amsoltimelimit
    )
    logger.debug("Running AMSOL7.1: SM5.42R (in hexadecane solvent)")
    run_external_command(
        f"{AMSOLEXE} < temp.in-hex > temp.o-hex", timeout=amsoltimelimit
    )

    # Extract data from AMSOL7.1 output files
    logger.debug("Extracting data from AMSOL7.1 water and hexadecane output files:")
    process_amsol_mol2(
        filenamewat=f"{current_path}/temp.o-wat",
        filenamehex=f"{current_path}/temp.o-hex",
        mol2file=temp_file,
        outputprefix=f"{current_path}/output",
    )

    # Clean
    for temp_file in ["fort.*", "temp.*", "*.log", "temp-working.mol2"]:
        for file in Path(current_path).glob(temp_file):
            os.remove(file)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        logger.debug("Usage: python calc_charge_solvation.py <mol2file>")
        sys.exit()
    calc_charge_solvation(sys.argv[1])
