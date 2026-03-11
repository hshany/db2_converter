import argparse
from pathlib import Path
from db2_converter.utils.conf_sample import available_methods

RMSthres = 0.5

def parse_args(version, header):
    parser = argparse.ArgumentParser(
        prog=f"db2_converter v{version}",
        description=header,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # IO group
    iogroup = parser.add_argument_group("io", "Options to define input and output")
    iogroup.add_argument(
        "-i",
        "--insmi",
        type=Path,
        default="input.smi",
        help='SMILES file with "SMILES LIG_NAME (conf_num)" on each line',
        required=True,
    )
    iogroup.add_argument(
        "--workingpath", type=Path, default=Path("."), help="Working path"
    )
    iogroup.add_argument(
        "--outputpath",
        type=Path,
        default=Path("."),
        help="Output path for mol2 and db2 files",
    )

    # Maxconfs group
    # The max_conf priority order:
    ## max_conf in smi file > max_conf in argument
    # The truncation order (with keep_max_conf, no truncation happen):
    ## rotateh -> limitconf
    maxconf_group = parser.add_argument_group("Maxconfs", "Options to control maxconfs")
    maxconf_group.add_argument(
        "-n", "--max_conf", type=int, default=600, help="Max number of conformers"
    )
    maxconf_group.add_argument(
        "--keep_max_conf",
        action="store_true",
        default=False,
        help="Do not alter max_conf according to settings in build3d.sh @ TLDR",
    )
    maxconf_group.add_argument(
        "-nr",
        "--max_conf_noring",
        type=int,
        default=30,
        help="Max conformers without rings, only valid without --keep_max_conf",
    )
    maxconf_group.add_argument(
        "--limitconf",
        action="store_true",
        default=False,
        help="Set limited conformational sampling, min(6*rot+3*ali_ring, max_conf)",
    )

    # Sample options group
    sampleopt_group = parser.add_argument_group(
        "SampleOptions", "Options for conformational sampling methods"
    )
    sampleopt_group.add_argument(
        "--sampletp",
        action="store_true",
        default=False,
        help="Determine protonation/tautomerization state",
    )
    sampleopt_group.add_argument(
        "-m",
        "--method",
        dest="samplopt",
        type=str,
        choices=available_methods,
        default="rdkit",
        help="Sampling method",
    )
    sampleopt_group.add_argument(
        "-c",
        "--checkstereo",
        action="store_true",
        help="Filter stereochemistry based on SMILES consistency",
    )
    sampleopt_group.add_argument(
        "-pb",
        "--PBfilter",
        action="store_true",
        help="Filter implausible conformers using PoseBusters",
    )
    sampleopt_group.add_argument(
        "-f",
        "--MMFFopt",
        action="store_true",
        default=False,
        help="Use MMFF to optimize conformers",
    )
    sampleopt_group.add_argument(
        "--cluster",
        action="store_true",
        help=f"Apply RMSD-based clustering (threshold: {RMSthres} Å) to conformers",
    )
    sampleopt_group.add_argument(
        "--bcl_option",
        type=str,
        default="",
        choices=["", "bcl-noBA", "bcl-noRing", "bcl-noBAnoRing"],
        help="BCL::Conf additional options",
    )
    sampleopt_group.add_argument(
        "--confgenx_option",
        type=str,
        default="",
        choices=["", "confgenx-OPLS2005", "confgenx-SOPLS", "confgenx-enableAmide"],
        help="ConfGenX additional options",
    )
    # sampleopt_group.add_argument(
    #     "-s",
    #     "--sampletp",
    #     action="store_true",
    #     default=False,
    #     help="sample tautomer and protomer, not supported now",
    # )

    # Rigid body group
    rigidbody_group = parser.add_argument_group(
        "RigidBody", "Options for configuration of rigid bodies"
    )
    rigidbody_group.add_argument(
        "--extra-fragsindex",
        nargs="+",
        type=int,
        default=[],
        help="One extra user-defined rigid fragment indices",
    )
    rigidbody_group.add_argument(
        "--extra-fragsmarts",
        type=str,
        default="",
        help="One extra user-defined fragment SMARTS",
    )
    # Only one of extra-fragsindex and extra-fragsmarts can be adopted (fragsindex > fragsmarts).
    rigidbody_group.add_argument(
        "--onlyextrafrags",
        action="store_true",
        default=False,
        help="Only use the extra fragments as rigid part",
    )
    rigidbody_group.add_argument(
        "--covalent",
        action="store_true",
        default=False,
        help="Use Si-centered rigid fragment (Si + heavy atoms within 2 bonds)",
    )

    # Mol2db2 group
    mol2db2_group = parser.add_argument_group("Mol2db2", "Options for mol2db2")
    mol2db2_group.add_argument(
        "--dock38",
        action="store_true",
        default=False,
        help="Generate DOCK3.8-formatted db2",
    )
    mol2db2_group.add_argument(
        "--reseth",
        action="store_true",
        default=False,
        help="Reset terminal polar hydrogens",
    )
    mol2db2_group.add_argument(
        "--rotateh",
        action="store_true",
        default=False,
        help="Rotate terminal polar hydrogens",
    )

    # Other options
    parser.add_argument(
        "--rerun",
        action="store_true",
        default=False,
        help="Remove previous running directories and rerun for the whole SMILES file",
    )
    parser.add_argument(
        "--mergeiso",
        action="store_true",
        default=False,
        help="Set the same name for different isomers in db2block output",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Set logger level at DEBUG (default is INFO))",
    )

    args = parser.parse_args()
    return args
