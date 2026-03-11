import os
import subprocess
import shutil
from pathlib import Path
import gzip
import logging

logger = logging.getLogger("db2_converter")

from db2_converter.pipeline import gen_conf
from db2_converter.utils.utils import check_type
from copy import deepcopy

def db2_converter(
    line,
    params,
    debug=True # haha
):
    workingpath = params["workingpath"]
    outputpath = params["outputpath"]
    max_conf = params["max_conf"]

    os.chdir(workingpath)
    lsp = line.split()
    zinc = lsp[1]
    if len(lsp) == 3 and check_type(lsp[2], int): # smi,name,max_conf
        max_conf = check_type(lsp[2], int)

    oneparams = deepcopy(params)
    oneparams["max_conf"] = max_conf
    for param in ["insmi","workingpath","outputpath","rerun","verbose"]:
        del oneparams[param]

    with open(f"{zinc}.smi", "w") as f:
        f.write(line + "\n")
    result_gen_conf = gen_conf(
        zinc=zinc,
        faillist=[],
        **oneparams,
    )
    os.chdir(workingpath)
    # collect result files (db2 and mol2) to outputpath
    if Path(f"{zinc}/all.db2.gz").exists():
        data = gzip.open(f"{zinc}/all.db2.gz", "rt").readlines()
        if data:
            shutil.move(f"{zinc}/all.db2.gz", f"{outputpath}/{zinc}.db2.gz")
            shutil.move(
                f"{zinc}/conformer.{zinc}.fixed.mol2",
                f"{outputpath}/conformer.{zinc}.fixed.mol2",
            )
    if not debug:
        subprocess.run(f"rm -r {zinc} {zinc}.smi", shell=True)

    return result_gen_conf

