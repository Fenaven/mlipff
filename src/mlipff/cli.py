import argparse


def add_convert_subparser(subparsers):
    """Add the 'convert' subparser."""
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
        "-i",
        dest="input_file",
        required=True,
        help="Input file path",
    )
    convert_parser.add_argument(
        "--output",
        "-o",
        dest="output_file",
        required=True,
        help="Output file path",
    )
    convert_parser.add_argument(
        "--replace_types",
        dest="replace_file",
        help="File for type replacement",
    )


def add_cut_data_subparser(subparsers):
    """Add the 'cut_data' subparser."""
    cut_data_parser = subparsers.add_parser(
        "cut_data",
        help="Cut data from system.data using nbh.dump",
    )
    cut_data_parser.add_argument(
        "--input",
        "-i",
        dest="input_file",
        required=True,
        help="Input data file path.",
    )
    cut_data_parser.add_argument(
        "--dump",
        "-d",
        dest="dump_file",
        required=True,
        help="Neighbourhood dump file path.",
    )
    cut_data_parser.add_argument(
        "--output",
        "-o",
        dest="output_file",
        required=True,
        help="Output data file path.",
    )


def add_cut_ff_subparser(subparsers):
    """Add the 'cut_ff' subparser."""
    cut_ff_parser = subparsers.add_parser(
        "cut_ff",
        help="Partition force field file",
    )
    cut_ff_parser.add_argument(
        "--input",
        "-i",
        dest="input_file",
        default="system.in.settings",
        help="Input file to read from.",
    )
    cut_ff_parser.add_argument(
        "--style",
        "-s",
        dest="style",
        default="lj/cut/coul/cut",
        help="Pair style to use in output.",
    )
    cut_ff_parser.add_argument(
        "--units",
        "-u",
        dest="units",
        default="metal",
        help="Units to use in output.",
    )


def add_create_input_subparser(subparsers):
    """Add the 'create_input' subparser."""
    input_parser = subparsers.add_parser(
        "create_input",
        help="Create Orca or LAMMPS input file using .xyz file and input template",
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
        "-i",
        dest="input_file",
        required=True,
        help="Input file path.",
    )
    input_parser.add_argument(
        "--xyz",
        dest="xyz_file",
        required=True,
        help="XYZ file path.",
    )
    input_parser.add_argument(
        "--output",
        "-o",
        dest="output_file",
        required=True,
        help="Output file path.",
    )


def add_cut_nbh_subparser(subparsers):
    """Add the 'cut_nbh' subparser."""
    cut_nbh_parser = subparsers.add_parser(
        "cut_nbh",
        help="Cut molecule using a sphere selection. Cutting from cfg or dump requires element_types file present",
    )
    cut_nbh_parser.add_argument(
        "--input",
        "-i",
        dest="input_file",
        required=True,
        help="Input file path, including extension.",
    )
    cut_nbh_parser.add_argument(
        "--output",
        "-o",
        dest="output_base",
        required=True,
        help="Output file base name.",
    )
    cut_nbh_parser.add_argument(
        "--coords",
        "--coordinates",
        dest="atom_coords",
        type=str,
        required=True,
        help="Comma-separated coordinates of the sphere center x,y,z.",
    )
    cut_nbh_parser.add_argument(
        "--radius",
        "-r",
        dest="radius",
        type=float,
        required=True,
        help="Radius of the selection sphere.",
    )
    cut_nbh_parser.add_argument(
        "--cut",
        dest="cut",
        action="store_true",
        help="If True, then cut molecules and add H to xyz file.",
    )
    cut_nbh_parser.add_argument(
        "--save-xyz",
        dest="save_xyz",
        action="store_true",
        help="Save extracted neighborhood in XYZ format.",
    )
    cut_nbh_parser.add_argument(
        "--save-pdb",
        dest="save_pdb",
        action="store_true",
        help="Save extracted neighborhood in PDB format.",
    )


def add_cut_dump_subparser(subparsers):
    """Add the 'cut_dump' subparser."""
    split_dump_parser = subparsers.add_parser(
        "cut_dump",
        help="Cut a dump file into separate files, one per timestep.",
    )
    split_dump_parser.add_argument(
        "--input",
        "-i",
        dest="input_file",
        required=True,
        help="Path to the input dump file.",
    )
    split_dump_parser.add_argument(
        "--output_prefix",
        "-o",
        dest="output_prefix",
        default="frame_",
        help="Prefix for the output files (default: 'frame_').",
    )


def add_cut_random_nbh_subparser(subparsers):
    """Add the 'cut_random_nbh' subparser."""
    cut_random_nbh_parser = subparsers.add_parser(
        "cut_random_nbh",
        help="Process MD frames by selecting a random atom and extracting a neighborhood.",
    )
    cut_random_nbh_parser.add_argument(
        "--frames",
        "-f",
        dest="num_frames",
        type=int,
        required=True,
        help="Number of frames to process.",
    )
    cut_random_nbh_parser.add_argument(
        "--prefix",
        "-p",
        dest="file_prefix",
        required=True,
        help="Prefix for the dump files (e.g., 'dump_').",
    )
    cut_random_nbh_parser.add_argument(
        "--radius",
        "-r",
        dest="radius",
        type=float,
        required=True,
        help="Radius for the neighborhood selection.",
    )
    cut_random_nbh_parser.add_argument(
        "--cut",
        "-c",
        dest="cut",
        action="store_true",
        help="Enable molecule cutting in the extracted neighborhood.",
    )


def add_filter_by_grade_subparser(subparsers):
    """Add the 'filter_by_grade' subparser."""
    filter_parser = subparsers.add_parser(
        "filter_by_grade",
        help="Filter configurations based on a minimum grade threshold.",
    )
    filter_parser.add_argument(
        "--input",
        "-i",
        dest="input_file",
        required=True,
        help="Input file path.",
    )
    filter_parser.add_argument(
        "--output",
        "-o",
        dest="output_file",
        required=True,
        help="Output file path.",
    )
    filter_parser.add_argument(
        "--threshold",
        "-t",
        dest="threshold",
        type=float,
        required=True,
        help="Minimum grade threshold.",
    )


def get_parser():
    """Create and return the argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        description="CLI for using MLIP-FF related functions"
    )
    subparsers = parser.add_subparsers(dest="command", help="Sub-command to run")

    # Register subcommands using functions
    subcommand_functions = [
        add_convert_subparser,
        add_cut_data_subparser,
        add_cut_ff_subparser,
        add_create_input_subparser,
        add_cut_nbh_subparser,
        add_cut_dump_subparser,
        add_cut_random_nbh_subparser,
        add_filter_by_grade_subparser,
    ]

    for func in subcommand_functions:
        func(subparsers)

    return parser


def parse_arguments():
    """Parse command-line arguments and return the parsed object."""
    return get_parser().parse_args()
