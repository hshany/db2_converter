#!/usr/bin/env python3

import os
import sys
import subprocess
import glob
import logging
import argparse
from pathlib import Path
import shutil
import time
from collections import defaultdict
from rdkit import Chem

def setup_logging():
    """Configure logging settings"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        stream=sys.stdout
    )

def run_command(cmd, working_dir, error_msg, stream_output=False):
    """Execute a shell command and handle errors"""
    try:
        if stream_output:
            # Use Popen for real-time output streaming
            process = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # Redirect stderr to stdout
                text=True,
                bufsize=1,  # Line buffered
                cwd=working_dir
            )

            # Stream output in real-time
            for line in process.stdout:
                print(line, end='')  # end='' because line already contains newline

            # Wait for process to complete and check return code
            return_code = process.wait()
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, cmd)

        else:
            # Use run for commands where we don't need to stream output
            logging.info(f"Running command: {cmd}")
            result = subprocess.run(
                cmd,
                check=True,
                shell=True,
                capture_output=True,
                text=True,
                cwd=working_dir
            )

        logging.info(f"Successfully ran command: {cmd}")

    except subprocess.CalledProcessError as e:
        logging.error(f"{error_msg}: {e}")
        if not stream_output and e.stderr:
            logging.error(f"Error output: {e.stderr}")
        raise

def convert_smi_to_mae(input_smi, output_mae, schrodinger_path, working_dir):
    """Convert .smi file to .mae format"""
    cmd = f"{schrodinger_path}/utilities/structconvert {input_smi} {output_mae}"
    run_command(cmd, working_dir, "Failed to convert .smi to .mae")
    return output_mae

def convert_mae_to_smi(input_mae, output_smi, schrodinger_path, working_dir):
    """Convert .mae file back to .smi format"""
    cmd = f"{schrodinger_path}/utilities/structconvert {input_mae} {output_smi}"
    run_command(cmd, working_dir, "Failed to convert .mae to .smi")
    return output_smi

def run_ligprep(input_smi, output_smi, schrodinger_path, working_dir, ph=7.4, pht=1.0):
    """Run Ligprep for protomer/tautomer/stereoisomer enumeration"""
    cmd = f"{schrodinger_path}/ligprep -epik -ph {ph} -pht {pht} -ismi {input_smi} -osmi {output_smi} -WAIT"
    run_command(cmd, working_dir, "Failed to run Epik")
    return output_smi

def run_build_ligand(
    input_smi, working_dir, n_conf=600, n_rot=30, covalent=False, cluster=False
):
    """Run build_ligand tool"""
    covalent_flag = " --covalent" if covalent else ""
    cluster_flag = " --cluster" if cluster else ""
    cmd = (f"OMP_NUM_THREADS=1 build_ligand -i {input_smi} -n {n_conf} -nr {n_rot} "
           f"--keep_max_conf -c -m confgenx{covalent_flag}{cluster_flag} --workingpath {working_dir} "
           f"--outputpath {working_dir}")
    run_command(cmd, working_dir, "Failed to run build_ligand", stream_output=True)

def concatenate_db2_files(working_dir, output_file):
    """Concatenate all db2.gz files into a single file"""
    db2_files = glob.glob(os.path.join(working_dir, "*.db2.gz"))
    if not db2_files:
        logging.warning("No db2.gz files found to concatenate")
        return
    cmd = f"cat {' '.join(db2_files)} > {output_file}"
    run_command(cmd, working_dir, "Failed to concatenate db2.gz files", stream_output=True)

def rename_smi(input_smi):
    name2smi = defaultdict(list)
    for line in open(input_smi).readlines():
        smi, name = line.split()
        name2smi[name].append(smi)
    with open(input_smi, 'w') as f:
        for name, smilist in name2smi.items():
            for i, smi in enumerate(smilist):
                f.write(f'{smi} {name}.{i}\n')

def filter_compounds(input_smi, output_smi):
    """
    Filter compounds from a .smi file, keeping only those with C, H, O, N, S, P and halogens.

    Parameters:
    input_smi (str): Path to input .smi file
    output_smi (str): Path to output .smi file

    Returns:
    tuple: (number of input compounds, number of filtered compounds)
    """
    # Define allowed elements (including halogens)
    allowed_elements = set(['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Si'])

    # Counter for statistics
    total_count = 0
    filtered_count = 0

    # Open output file for writing
    with open(output_smi, 'w') as outfile:
        # Read input file line by line
        with open(input_smi, 'r') as infile:
            for line in infile:
                total_count += 1

                # Split line into SMILES and name (if present)
                parts = line.strip().split()
                smiles = parts[0]
                name = parts[1] if len(parts) > 1 else ""

                # Try to create molecule from SMILES
                mol = Chem.MolFromSmiles(smiles)

                if mol is None:
                    continue

                # Get all unique elements in the molecule
                atoms = set(atom.GetSymbol() for atom in mol.GetAtoms())

                # Check if all elements are in allowed set
                if atoms.issubset(allowed_elements):
                    filtered_count += 1
                    # Write to output file, maintaining the same format as input
                    if name:
                        outfile.write(f"{smiles}\t{name}\n")
                    else:
                        outfile.write(f"{smiles}\n")

    return total_count, filtered_count

def main():
    parser = argparse.ArgumentParser(description="Process .smi files for docking")
    parser.add_argument("input_smi", help="Input .smi file")
    parser.add_argument("--schrodinger", help="Path to SCHRODINGER directory", 
                       default="/home/hehuang/soft_shared/schrodinger2023-2")
    parser.add_argument("--working_path", help="Working directory path", 
                       default=os.getcwd())
    parser.add_argument("--output_db2", help="Output db2.gz file name",
                       default="combined_output.db2.gz")
    parser.add_argument("--ph", type=float, default=7.4, 
                       help="pH value for Epik (default: 7.4)")
    parser.add_argument("--pht", type=float, default=1.0, 
                       help="pH tolerance for Epik (default: 1.0)")
    parser.add_argument("--n_conf", type=int, default=600, 
                       help="Number of conformers (default: 600)")
    parser.add_argument("--n_rot", type=int, default=30, 
                       help="Number of rotatable bonds (default: 30)")
    parser.add_argument("--covalent", action="store_true",
                       help="Enable covalent mode (SiH3-based alignment)")

    args = parser.parse_args()

    # Setup logging
    setup_logging()

    # Create working directory if it doesn't exist
    working_path = Path(args.working_path)
    working_path.mkdir(parents=True, exist_ok=True)

    try:
        # Handle input file and working directory
        input_smi_path = Path(args.input_smi).absolute()
        working_path = Path(args.working_path).absolute()
        input_smi_name = input_smi_path.name
        working_smi = working_path / input_smi_name

        # Only copy if input file is not already in working directory
        if input_smi_path != working_smi:
            if working_smi.exists():
                logging.warning(f"File {working_smi} already exists in working directory")
            else:
                shutil.copy2(input_smi_path, working_smi)

        # Set up intermediate file paths in working directory
        ligprep_smi = Path(input_smi_name.replace('.smi', '.prep.smi'))
        filtered_smi = Path(ligprep_smi.name.replace('.smi', '.filtered.smi'))

        # Run Ligprep
        run_ligprep(input_smi_name, ligprep_smi, args.schrodinger, str(working_path), args.ph, args.pht)

        # Filter out problematic compounds
        filter_compounds(ligprep_smi, filtered_smi)

        if os.path.exists(filtered_smi):
            # Rename protomer/tautomer
            rename_smi(filtered_smi)

            # Run build_ligand
            run_build_ligand(
                str(filtered_smi),
                str(working_path),
                args.n_conf,
                args.n_rot,
                args.covalent,
                cluster=True,
            )

            # Concatenate output files
            output_db2 = args.output_db2
            concatenate_db2_files(str(working_path), output_db2)

            logging.info(f"Successfully created combined db2.gz file: {output_db2}")
        else:
            logging.info(f"No valid protomers generated by Ligprep")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()
