import shutil
import logging
from datetime import datetime
from tabulate import tabulate
from db2_converter.utils.utils import exist_size, run_external_command

logger = logging.getLogger("config")
config = {
    "all":{
        "BABEL_EXE": f"obabel",
        "AMSOLEXE": f"{__file__.split('config.py')[0]}/amsol/amsol7.1_patched",
        "ANTECHAMBER": f"antechamber",
        "UNICON_EXE": f"/home/hehuang/soft_shared/unicon/unicon",
    },
    "conformator":{
        "CONF_EXE": f"conformator"
    },
    "bcl":{
        "BCLBASE": f"/your/path/to/bcl-4.2.0-Linux-x86_64",
        "BCL": f"/your/path/to/bcl-4.2.0-Linux-x86_64/bcl.exe molecule:ConformerGenerator"
    },
    "confgenx":{
        "CONFGENX": "/home/hehuang/soft_shared/schrodinger2023-2/confgenx",
        "SCHUTILS": "/home/hehuang/soft_shared/schrodinger2023-2/utilities"
    },
    "ccdc":{
        "CCDC_PYTHON3": f"/your/path/to/ccdc-software/csd-python-api/miniconda/bin/python3.9",
        "CCDC_pyscript": f"{__file__.split('config.py')[0]}/utils/ccdc_confgen.py",
        "CCDC_activate_code": "xxx", # only can be used privately!!!
        "CCDC_activate_command": f"/your/path/to/ccdc-utilities/software-activation/bin/ccdc_activator -a ",
        "CCDC_activate_node":"xxx"
    }
}


def check_config_avail(samplopt):
    header = ["level","name","value","available"]
    lines = []
    methods = samplopt.split("_")
    exclude_keys = ["BCL","CCDC_activate_code","CCDC_activate_command","CCDC_activate_node"]
    for key in ["all"] + methods:
        if key in config:
            for subkey in config[key]:
                avail = True
                value = config[key][subkey]
                if subkey in exclude_keys:
                    avail = "/"
                else:
                    if not exist_size(config[key][subkey]):
                        if shutil.which(value):
                            value = shutil.which(value)
                        else:
                            avail = False
                lines.append([key,subkey,value,avail])
    logger.info(
        "Config availability:\n"
        + tabulate(tabular_data=lines, headers=header,tablefmt="grid")
        + "\n"
    )
    not_avail_lines = [ line for line in lines if line[-1] == False]
    if not_avail_lines:
        logger.error("Not all configs are set up!!!\n")
        return
    return True


def check_UNICON_license():
    UNICON_EXE = config["all"]["UNICON_EXE"]
    result = run_external_command(UNICON_EXE)
    result_stdout_lines = result.stdout.split("\n")
    license_line = [ line for line in result_stdout_lines if line.startswith("Your license is valid until:")][0]
    date_str = license_line.split("Your license is valid until:")[1].strip()
    if date_str:
        date_format = "%A, %d %B %Y"
        license_expiry_date = datetime.strptime(date_str, date_format)
        current_date = datetime.now()
        if license_expiry_date > current_date:
            return True
        else:
            logger.error("!!! UNICON license expired. !!!")
    else:
        logger.error("!!! UNICON license expired. !!!")
