from utils import (
    convert,
    cut_data,
    cut_ff,
    cut_nbh,
    generate_orca_input,
    modify_lammps_data,
    cut_dump,
)
from cli import parse_arguments, get_parser


class UnknownCommandError(Exception):
    """Exception raised when an unknown CLI command is used."""

    pass


def handle_convert(args):
    convert(
        args.from_format,
        args.to_format,
        args.input_file,
        args.output_file,
        args.replace_file,
    )


def handle_cut_data(args):
    cut_data(args.input_file, args.dump_file, args.output_file)


def handle_cut_ff(args):
    cut_ff(args.input_file, args.style, args.units)


def handle_cut_nbh(args):
    atom_coords = list(map(float, args.atom_coords.split(",")))
    cut_nbh(args.input_file, args.output_base, atom_coords, args.radius, args.cut)


def handle_create_input(args):
    if args.orca:
        generate_orca_input(args.input_file, args.xyz_file)
    elif args.lmp:
        modify_lammps_data(args.input_file, args.xyz_file)


def handle_cut_dump(args):
    cut_dump(args.input_file, args.output_prefix)


# Mapping CLI commands to handler functions
COMMAND_HANDLERS = {
    "convert": handle_convert,
    "cut_data": handle_cut_data,
    "cut_ff": handle_cut_ff,
    "cut_nbh": handle_cut_nbh,
    "create_input": handle_create_input,
    "cut_dump": handle_cut_dump,
}


def main():
    args = parse_arguments()

    if args.command in COMMAND_HANDLERS:
        COMMAND_HANDLERS[args.command](args)
    else:
        raise UnknownCommandError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    try:
        main()
    except UnknownCommandError as e:
        print(f"Error: {e}")
        get_parser().print_help()
        exit(1)
