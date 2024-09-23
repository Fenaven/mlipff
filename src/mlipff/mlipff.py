from utils import convert, cut_data, cut_ff, cut_nbh, generate_orca_input
from cli import parse_arguments, get_parser

def main():
    args = parse_arguments()

    if args.command == "convert":
        convert(
            args.from_format,
            args.to_format,
            args.input_file,
            args.output_file,
        )
    elif args.command == "cut_data":
        cut_data(
            args.data_file,
            args.dump_file,
            args.output_file,
        )
    elif args.command == "cut_ff":
        cut_ff(
            args.input_file,
            args.style,
        )
    elif args.command == "cut_nbh":
        atom_coords = list(map(float, args.atom_coords.split(",")))
        cut_nbh(
            args.input_file,
            args.output_base,
            atom_coords,
            args.radius,
            args.cut
        )
    elif args.command == "create_input":
        generate_orca_input(args.xyz_file)
    else:
        # If no command is provided, print the help message
        get_parser().print_help()

if __name__ == "__main__":
    main()
