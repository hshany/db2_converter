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


def convert_sdf_to_mol2(insdf, outmol2, tool="auto"):
    """
    Convert SDF to MOL2 using specified tool.

    Args:
        insdf: Input SDF file path
        outmol2: Output MOL2 file path
        tool: Conversion tool to use. Options:
            "auto" (default) - Try openbabel first, fallback to schrodinger if fails
            "openbabel" - Use OpenBabel (free, recommended)
            "schrodinger" - Use Schrodinger structconvert (requires license)

    Raises:
        RuntimeError: If conversion fails or no suitable tool available
    """
    import logging
    import shutil

    logger = logging.getLogger("db2_converter")

    # Allow config override of default tool
    if tool == "auto":
        tool = config.get("all", {}).get("SDF2MOL2_TOOL", "auto")

    # Try OpenBabel
    if tool in ("auto", "openbabel"):
        obabel = shutil.which(config.get("all", {}).get("BABEL_EXE", "obabel"))
        if obabel:
            logger.debug(f"Converting {insdf} -> {outmol2} using OpenBabel")
            result = run_external_command(
                f"{obabel} -isdf {insdf} -omol2 -O {outmol2}",
                stderr=subprocess.STDOUT
            )
            if result.returncode == 0 and exist_size(outmol2):
                return
            else:
                logger.warning(f"OpenBabel conversion failed (exit code {result.returncode})")
                if tool == "openbabel":
                    raise RuntimeError(f"OpenBabel conversion failed for {insdf}")
                # If auto mode, try schrodinger fallback
        else:
            logger.warning("OpenBabel (obabel) not found in PATH")
            if tool == "openbabel":
                raise RuntimeError(
                    "OpenBabel not available. Install via: conda install -c conda-forge openbabel"
                )

    # Try Schrodinger structconvert
    if tool in ("auto", "schrodinger"):
        if "confgenx" in config and "SCHUTILS" in config["confgenx"]:
            structconvert_path = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")
            if exist_size(structconvert_path) or shutil.which("structconvert"):
                structconvert = structconvert_path if exist_size(structconvert_path) else "structconvert"
                logger.debug(f"Converting {insdf} -> {outmol2} using Schrodinger structconvert")
                result = run_external_command(
                    f"{structconvert} {insdf} {outmol2}",
                    stderr=subprocess.STDOUT
                )
                if result.returncode == 0 and exist_size(outmol2):
                    return
                else:
                    logger.warning(f"structconvert conversion failed (exit code {result.returncode})")
                    raise RuntimeError(f"Schrodinger structconvert failed for {insdf}")
            else:
                logger.warning("Schrodinger structconvert not found")

        if tool == "schrodinger":
            raise RuntimeError(
                "Schrodinger structconvert not available. "
                "Configure SCHUTILS path in config or use tool='openbabel'"
            )

    # No tool succeeded
    raise RuntimeError(
        f"SDF to MOL2 conversion failed for {insdf}.\n"
        f"Attempted tool: {tool}\n"
        f"Please install OpenBabel: conda install -c conda-forge openbabel\n"
        f"Or configure Schrodinger tools in config"
    )
