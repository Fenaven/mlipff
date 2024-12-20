import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description="Tool for converting and manipulating configuration files."
    )
    subparsers = parser.add_subparsers(
        dest="command",
        help="Sub-command to run",
    )

    # Convert command
    convert_parser = subparsers.add_parser(
        "convert",
        help="Convert files between formats",
    )
    convert_parser.add_argument(
        "--from",
        dest="from_format",
        choices=["dump", "orca", "orcadump"],
        required=True,
        help="Input format",
    )
    convert_parser.add_argument(
        "--to",
        dest="to_format",
        choices=["cfg"],
        required=True,
        help="Output format",
    )
    convert_parser.add_argument(
        "--input",
        dest="input_file",
        required=True,
        help="Input file path",
    )
    convert_parser.add_argument(
        "--output",
        dest="output_file",
        required=True,
        help="Output file path",
    )

    convert_parser.add_argument(
        "--replace_types",
        dest="replace_file",
        help="File for type replacement",
    )

    # Cut_data command
    cut_data_parser = subparsers.add_parser(
        "cut_data",
        help="Cut data from system.data using nbh.dump",
    )
    cut_data_parser.add_argument(
        "--input",
        dest="input_file",
        required=True,
        help="Input data file path.",
    )
    cut_data_parser.add_argument(
        "--nbh",
        dest="dump_file",
        required=True,
        help="Neighbourhood dump file path.",
    )
    cut_data_parser.add_argument(
        "--output",
        dest="output_file",
        required=True,
        help="Output data file path.",
    )

    # Cut_ff command
    cut_ff_parser = subparsers.add_parser(
        "cut_ff",
        help="Partition force field file",
    )
    cut_ff_parser.add_argument(
        "--input",
        dest="input_file",
        default="system.in.settings",
        help="Input file to read from.",
    )
    cut_ff_parser.add_argument(
        "--style",
        dest="style",
        default="lj/cut/coul/cut",
        help="Pair style to use in output.",
    )
    cut_ff_parser.add_argument(
        "--units",
        dest="units",
        default="metal",
        help="Units to use in output.",
    )

    # Create_input command
    input_parser = subparsers.add_parser(
        "create_input",
        help="Create Orca input file using .xyz file",
    )
    input_parser.add_argument(
        "--orca",
        action="store_true",
        help="Create input for Orca",
    )
    input_parser.add_argument(
        "--lmp",
        action="store_true",
        help="Create input for LAMMPS",
    )
    input_parser.add_argument(
        "--input",
        dest="input_file",
        required=True,
        help="input file path.",
    )
    input_parser.add_argument(
        "--xyz",
        dest="xyz_file",
        required=True,
        help="xyz file path.",
    )

    # Cut_nbh command
    cut_nbh_parser = subparsers.add_parser(
        "cut_nbh", help="Cut molecule using a sphere selection"
    )
    cut_nbh_parser.add_argument(
        "--input",
        dest="input_file",
        required=True,
        help="Input dump file path.",
    )
    cut_nbh_parser.add_argument(
        "--output",
        dest="output_base",
        required=True,
        help="Output file base name.",
    )
    cut_nbh_parser.add_argument(
        "--coords",
        dest="atom_coords",
        required=True,
        help="Comma-separated coordinates of the sphere center (x,y,z).",
    )
    cut_nbh_parser.add_argument(
        "--radius",
        dest="radius",
        type=float,
        required=True,
        help="Radius of the selection sphere.",
    )

    cut_nbh_parser.add_argument(
        "--cut",
        action="store_true",
        help="If True, then cut molecules and add H to xyz file.",
    )

    return parser


def parse_arguments():
    return get_parser().parse_args()
