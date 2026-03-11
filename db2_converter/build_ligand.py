db2_converter_version = "0.3"

import os
import shutil
from pathlib import Path
import logging
from functools import partial
import time
from tabulate import tabulate
import gzip

# disable rdkit warning log
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.warning")

from db2_converter.parse_args import parse_args
from db2_converter.pipeline import sample_tp_unicon, write_enumerated_smifile
from db2_converter.utils.utils import exist_size, next_mol2_lines
from db2_converter.mol2db2 import mol2db2
from db2_converter.utils.match_frags import mol_to_ring_frags
from db2_converter.utils.ccdc_confgen import activate_ccdc
from db2_converter.db2_converter import db2_converter
from db2_converter.config import check_config_avail, check_UNICON_license
from db2_converter.pipeline import sb_smarts

def create_header():
    from db2_converter import config as cfg

    header = ""
    header += "\n"
    header += "#" * 73
    header += r"""
  ____   ____  ____                                        _              
 |  _ \ | __ )|___ \    ___  ___   _ __ __   __ ___  _ __ | |_  ___  _ __ 
 | | | ||  _ \  __) |  / __|/ _ \ | '_ \\ \ / // _ \| '__|| __|/ _ \| '__|
 | |_| || |_) |/ __/  | (__| (_) || | | |\ V /|  __/| |   | |_|  __/| |   
 |____/ |____/|_____|  \___|\___/ |_| |_| \_/  \___||_|    \__|\___||_|   
"""
    header += "\n"
    header += f"Packaged from:\n {__file__}\n"
    header += f"Remember to modify the config file for the third-party softwares:\n {cfg.__file__}\n\n"
    header += "#" * 73 + "\n"
    return header


def main():
    header = create_header()
    args = parse_args(version=db2_converter_version, header=header)
    args.outputpath = args.outputpath.absolute()
    args.workingpath = args.workingpath.absolute()
    infile = args.insmi

    if args.rerun:
        for exist_dir in [args.workingpath, args.outputpath]:
            shutil.rmtree(exist_dir, ignore_errors=True)
    for dirname in [args.outputpath, args.workingpath]:
        dirname.mkdir(exist_ok=True)
    try:
        shutil.copy(infile, args.workingpath)
    except shutil.SameFileError:
        pass

    os.chdir(args.workingpath)

    # logging
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logmode = "w" if args.rerun else "a"
    logfile = args.workingpath / f"{infile.name}.{args.samplopt}.log"
    logging.basicConfig(
        level=log_level,
        format="%(levelname)7s:%(name)25s:[%(asctime).19s] %(message)s",
        handlers=[
            logging.FileHandler(logfile, logmode),
            logging.StreamHandler(),
        ],
    )
    logger = logging.getLogger("build_ligand")
    logger.info(header)
    logger.info(str(vars(args)) + "\n")

    # check parameters availability
    if not check_config_avail(args.samplopt):
        return
    if args.sampletp:
        if not check_UNICON_license():
            return
    # activate ccdc on targeted machine
    if "ccdc" in args.samplopt:
        activate_ccdc()

    params = vars(args)
    # print(params.keys())
    partial_db2_converter = partial(
        db2_converter,
        params=params,
    )

    if args.sampletp:
        infile = sample_tp_unicon(infile)
    lines = [line for line in infile.read_text().split("\n") if line]


    def oldrun(line,args): # To obtain information from generated files
        max_conf = args.max_conf
        lsp = line.split()
        if len(lsp) > 2:
            smi,name,max_conf = line.split()
            max_conf = int(max_conf)
        else:
            smi,name = line.split()
        gen_db2gzfile = args.outputpath / f"{name}.db2.gz"
        gen_mol2file = args.outputpath / f"conformer.{name}.fixed.mol2"
        if (exist_size(gen_db2gzfile) and exist_size(gen_mol2file)):
            telapsed = "/"
            ringmol = Chem.MolFromSmiles(smi)
            N_ring = len(mol_to_ring_frags(ringmol, cutsmarts=sb_smarts))
            NrotHs, multiplier = mol2db2.mol2db2_to_numhyds(
                smifile="", mol2file=str(gen_mol2file), removemol2=False
            )
            if not args.keep_max_conf and N_ring:
                if NrotHs in [4, 5]:
                    max_conf = max_conf // 30
                if NrotHs in [2, 3]:
                    max_conf = max_conf // 3
            N_max_conf_in = max_conf
            N_act_conf = len(list(next_mol2_lines(gen_mol2file)))
            N_act_conf_out_rotH = N_act_conf * multiplier
            with gzip.open(gen_db2gzfile, "rt") as f:
                N_db2_part = sum(1 for line in f if line.startswith("E"))
            return [telapsed,N_max_conf_in,N_act_conf,N_act_conf_out_rotH,N_ring,N_db2_part,""]
        else:
            return ["/"]*6 + [[]]


    def newrun(line): # run db2_converter workflow
        start = time.time()
        (
            faillist,
            N_max_conf_in,
            N_act_conf,
            N_act_conf_out_rotH,
            N_ring,
            N_db2_part,
        ) = partial_db2_converter(line=line)
        end = time.time()
        telapsed = "{:.4f}".format(end - start)
        return [telapsed,N_max_conf_in,N_act_conf,N_act_conf_out_rotH,N_ring,N_db2_part,faillist]

    allfaillist = []
    if args.checkstereo:
        enu_infile = f"{args.workingpath}/{infile.stem}.enumerated.smi"
        faillist = write_enumerated_smifile(
            lines, enu_infile, args.samplopt
        )
        if faillist: allfaillist = faillist
        lines = [line for line in Path(enu_infile).read_text().split("\n") if line]
    else:
        logger.info(
            ">>> Stereochemistry Enumeration is skipped. Be careful on your results!!! "
        )
        lines = [line for line in infile.read_text().split("\n") if line]
    names_done = []
    # Exclude generated files in this run
    lastname = ""
    if not args.rerun:
        if exist_size(logfile) and len(list((args.workingpath).glob("*.db2.gz"))):
            with open(logfile) as f:
                for line in f:
                    line = line[56:]
                    if line.startswith("############### Now dealing with"):
                        lastname = line.split("############### Now dealing with")[1].split("...")[0].strip()
                    if line.startswith(f"############### Finished with {lastname}"): # have finished
                        names_done.append(lastname)
                        lastname = ""
            print("names_done",names_done)
            lines = [ line for line in lines if line.split()[1] not in names_done ]
            print("lines",lines)
            if not lines:
                logger.info(">>> All DB2 generation task has completed in the last run. Exit...")
                return
            logger.info(">>> Continue last run...\n")

    allsum = []
    runleft = True if args.rerun else False
    for line in lines:
        lsp = line.split()
        smi = lsp[0]
        name = lsp[1]
        results = newrun(line)
        telapsed,N_max_conf_in,N_act_conf,N_act_conf_out_rotH,N_ring,N_db2_part,faillist = results
        failinfo = [onefail[-1] for onefail in faillist if onefail]
        if failinfo: failinfo = "_".join(failinfo)
        else: failinfo = ""
        allsum.append([
            name,
            args.samplopt,
            telapsed,
            N_max_conf_in,
            N_act_conf,
            N_act_conf_out_rotH,
            N_ring,
            N_db2_part,
            failinfo,
            smi,
        ])
        allfaillist.append(faillist)

    # Tabulate output
    table_header = [
        "Name",
        "method",
        "telapsed",
        "N_in",
        "N_out",
        "N_out_rotH",
        "N_ring",
        "N_db2part",
        "FAILED",
        "SMILES",
    ]
    logger.info(
        "ALL SUMMARY:\n"
        + tabulate(tabular_data=allsum, headers=table_header, tablefmt="grid")
    )

    # Write fail summary for a list of smis
    if allfaillist:
        writefaillist = ["\t".join(fail) for faillist in allfaillist if faillist for fail in faillist if fail]
        logger.info("<<<<< Writing faillist... >>>>>")
        with open(f"{args.outputpath}/{infile.name}.{args.samplopt}.faillist", "w") as f:
            f.write("\n".join(writefaillist))
            f.write("\n")


if __name__ == "__main__":
    main()
