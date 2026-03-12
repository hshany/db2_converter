import os
import shutil
import subprocess

from db2_converter.config import config
from db2_converter.utils.utils import (
    ATOMTYPE,
    exist_size,
    next_mol2_lines,
    run_external_command,
)


def swap_dummy_si_in_sdf(insdf, outsdf):
    si_atom_ids = []
    with open(insdf, "r") as f:
        lines = f.readlines()
    if len(lines) < 4:
        shutil.copy(insdf, outsdf)
        return si_atom_ids
    counts_line = lines[3]
    try:
        natoms = int(counts_line[0:3])
    except ValueError:
        natoms = int(counts_line.split()[0])
    start = 4
    end = min(start + natoms, len(lines))
    for idx in range(start, end):
        line = lines[idx]
        atom_id = idx - 3  # 1-based atom index
        if len(line) >= 34:
            elem = line[31:34].strip()
            if elem == "Si":
                si_atom_ids.append(atom_id)
                lines[idx] = f"{line[:31]} C {line[34:]}"
        else:
            tokens = line.split()
            if len(tokens) >= 4 and tokens[3] == "Si":
                si_atom_ids.append(atom_id)
                tokens[3] = "C"
                lines[idx] = " ".join(tokens) + "\n"
    with open(outsdf, "w") as f:
        f.writelines(lines)
    return si_atom_ids


def restore_dummy_si_in_mol2(mol2file, si_atoms):
    if not si_atoms:
        return
    mol2blocks = [x for x in next_mol2_lines(mol2file)]
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
            if atom_id in si_atoms:
                new_atompart.append(
                    ATOMTYPE.format(
                        atom_id,
                        si_atoms[atom_id],
                        float(items[2]),
                        float(items[3]),
                        float(items[4]),
                        "Si",
                        items[6],
                        items[7],
                        float(items[8]),
                    )
                )
            else:
                new_atompart.append(line)
        new_blocks.append("".join(startpart + new_atompart + bondpart))
    with open(mol2file, "w") as f:
        f.write("".join(new_blocks))


def convert_sdf_to_mol2(insdf, outmol2):
    structconvert = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")
    run_external_command(
        f"{structconvert} {insdf} {outmol2}",
        stderr=subprocess.STDOUT,
    )
