import subprocess
import signal
import os
from rdkit import Chem
import logging
logger = logging.getLogger("DB2 generation")

# Antechamber mol2 format
ATOMTYPE = (
    "    {:>3d} {:<8s} {:>10.4f} {:>10.4f} {:>10.4f} {:<6s}{:>6s} {:<6s}{:>12.6f}\n"
)
BONDTYPE = "{:>6d}{:>6d}{:>6d} {:<4s}\n"


def next_mol2_lines(infile):
    """Method to return one mol2 block once."""
    lines = list()
    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line:
            if len(lines) == 0:
                lines.append(line)
            else:  # in case there are multiple mol2blocks in infile
                yield lines
                lines = list()
                lines.append(line)
        else:
            lines.append(line)
    yield lines


def derive_first_mol2(inmol2, outmol2):
    with open(outmol2, "w") as ofp:
        start = True
        for line in open(inmol2):
            if "@<TRIPOS>MOLECULE" in line and start:
                start = False
                ofp.write(line)
            elif "@<TRIPOS>MOLECULE" in line and not start:
                break
            else:
                ofp.write(line)


def update_mol2block_from_mol(mol2block, newmol):
    p1 = mol2block.index("@<TRIPOS>ATOM\n")
    p2 = mol2block.index("@<TRIPOS>BOND\n")
    Natom, Nbond = mol2block[2].strip().split()[0:2]

    startpart = mol2block[:p1]
    atompart = ["@<TRIPOS>ATOM\n"]
    bondpart = ["@<TRIPOS>BOND\n"]
    new_pos = newmol.GetConformer().GetPositions()

    for cout, i in enumerate(range(p1 + 1, p1 + 1 + int(Natom))):
        items = mol2block[i].strip().split()
        x, y, z = new_pos[cout]
        atompart.append(
            ATOMTYPE.format(
                int(items[0]),
                items[1],
                float(x),
                float(y),
                float(z),
                items[5],
                items[6],
                items[7],
                float(items[-1]),
            )
        )

    for i in range(p2 + 1, p2 + 1 + int(Nbond)):
        items = mol2block[i].strip().split()
        bondpart.append(
            BONDTYPE.format(int(items[0]), int(items[1]), int(items[2]), items[3])
        )

    return startpart + atompart + bondpart + ["\n"]


def run_external_command(command_str,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        timeout=3600,
                        log=None):
    """
    Run an external command with a timeout, ensuring the entire process group
    is killed if any exception occurs during execution or if the process fails.
    """
    # Fallback to a default logger if none provided
    if log is None:
        import logging
        log = logging.getLogger(__name__)

    proc = None
    try:
        # We create a new session so we can kill the entire process group if needed
        proc = subprocess.Popen(
            command_str,
            stdout=stdout,
            stderr=stderr,
            shell=True,
            universal_newlines=True,
            start_new_session=True
        )

        out, err = proc.communicate(timeout=timeout)
        result = subprocess.CompletedProcess(
            args=command_str,
            returncode=proc.returncode,
            stdout=out,
            stderr=err
        )
        if result.stdout:
            log.debug(result.stdout)
        if result.stderr:
            log.error(result.stderr)

        return result

    except subprocess.TimeoutExpired as e:
        # If the process times out, log an error
        log.error(f"Command '{command_str}' timed out after {timeout} seconds.")
        raise

    except Exception as e:
        # Handle any other exceptions
        log.error(f"Command '{command_str}' failed with error: {str(e)}")
        raise

    finally:
        # Kill the process group only if proc exists and is still running
        if proc is not None and proc.poll() is None:
            try:
                log.info(f'Killing process group {proc.pid}')
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except Exception as kill_error:
                log.error(f"Failed to kill process group: {kill_error}")


def exist_size(testfile):
    return os.path.exists(testfile) and os.path.getsize(testfile)


def check_mol2_smi(mol2file, smi, isomericSmiles=False):
    try:
        probesmi = Chem.MolToSmiles(
            Chem.MolFromMol2File(mol2file), isomericSmiles=isomericSmiles
        )
        return probesmi == Chem.MolToSmiles(
            Chem.MolFromSmiles(smi), isomericSmiles=isomericSmiles
        )
    except:
        return False

def check_type(name, type):
    try:
        return type(name)
    except:
        return

def raise_errlog(error,logger):
    error_dict = {
        "1generate":
            ">>> No conformers sampled, please check your input smi.",
        "2fixmol2":
            ">>> No viable confomrers after fixmol2.",
        "3chemistrycheck":
            ">>> No viable conformers after chemistry check.",
        "4PBfilter":
            ">>> No viable conformers after PoseBusters filter.",
        "5amsolfail":
            ">>> AMSOL calculation of desolvation/charge failed...",
        "9nulldb2gz":
            ">>> Null db2gz! Null db2gz!! Null db2gz!!!"
    }
    logger.error(error_dict[error])
    logger.error("The task has failed! Generation stopped.\n")
