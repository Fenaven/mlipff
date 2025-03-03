import os
from typing import Set
from mlipff_configuration import Configuration
from ocellar import molecule


def convert(
    from_format: str,
    to_format: str,
    input_file: str,
    output_file: str,
    replace_file: str,
) -> None:
    """
    Convert a file from one format to another using the provided configuration.

    Parameters
    ----------
    from_format : str
        The input file format.
    to_format : str
        The desired output file format.
    input_file : str
        Path to the input file.
    output_file : str
        Path to the output file.
    replace_file : str
        Path to the replace file for configuration.
    """
    config = Configuration()

    if from_format == "dump":
        config.read_dump(input_file, True, replace_file)
    elif from_format == "orcadump":
        config.read_orca(input_file, True, replace_file)
    elif from_format == "orca":
        config.read_orca(input_file, False, replace_file)
    else:
        print(f"Unsupported input format: {from_format}")
        return

    if to_format == "cfg":
        config.write_cfg(output_file)
    else:
        print(f"Unsupported output format: {to_format}")
        return


def cut_data(input_file: str, dump_file: str, output_file: str) -> None:
    # ДОБАВИТЬ ВОЗМОЖНОСТЬ ЧИТАТЬ НЕ ИЗ ДАМПА, А ИЗ КФГ
    """
    Cut data from the data file based on atom IDs present in the dump file and write to the output file.

    Parameters
    ----------
    input_file : str
        Path to the data file.
    dump_file : str
        Path to the dump file.
    output_file : str
        Path to the output file.
    """

    def handle_section(
        line: str, section: str, tmp_content: list, ids: Set[int], atom_coords_map: dict
    ) -> None:
        """Handle the parsing of a section in the data file and update the content accordingly."""
        if (line.strip() == "") or (section == "header"):
            tmp_content.append(line.rstrip() + "\n")
            return

        splitted_line = line.split()
        type_set = type_sets[section]

        if section == "atoms":
            atom_id = int(splitted_line[0])
            if atom_id in ids:
                # Replace coordinates by those from dump
                splitted_line[-3:] = map(str, atom_coords_map[atom_id])
                tmp_content.append(" ".join(splitted_line) + "\n")
                counts[section] += 1
                type_set.add(int(splitted_line[2]))

        else:
            line_ids = list(map(int, splitted_line[2:]))
            if all(atom_id in ids for atom_id in line_ids):
                counts[section] += 1
                type_set.add(int(splitted_line[1]))
                tmp_content.append(line.rstrip() + "\n")

    ids: Set[int] = set()
    atom_coords_map = dict()
    type_sets = {
        t: set() for t in ["atoms", "bonds", "angles", "dihedrals", "impropers"]
    }
    counts = {
        t: 0 for t in type_sets
    }  # dict with keys atoms, bonds, ... for counting ids
    tmp_content = []
    current_section = "header"

    # Parse dump to collect info about ids and coordinates
    with open(dump_file, "r", encoding="utf-8") as dumpfile:
        for line in dumpfile:
            if "ITEM: NUMBER" in line:
                atom_number = int(next(dumpfile).strip())
            elif "ITEM: ATOMS" in line:
                coordinates_start_index = line.split()[2:].index("x")
                for _ in range(atom_number):
                    splitted_line = next(dumpfile).split()
                    pos_x, pos_y, pos_z = map(
                        float,
                        splitted_line[
                            coordinates_start_index : coordinates_start_index + 3
                        ],
                    )
                    ids.add(int(splitted_line[0]))
                    atom_coords_map[int(splitted_line[0])] = (pos_x, pos_y, pos_z)

    with (
        open(input_file, "r", encoding="utf-8") as datafile,
        open(output_file, "w", encoding="utf-8", newline="\n") as output,
    ):
        # Parse datafile and write directly to output
        for line in datafile:
            if not line.strip():
                tmp_content.append("\n")
                continue
            first_word = line.split()[0].lower()
            if first_word in counts:
                current_section = first_word
                tmp_content.append(line.rstrip() + "\n")
            else:
                handle_section(line, current_section, tmp_content, ids, atom_coords_map)

        # Update header information and write to output
        for i, line in enumerate(tmp_content):
            for name in counts:
                if name in line:
                    tmp_content[i] = f"     {counts[name]}  {name}\n"
                    break
            else:  # if there was no match in counts
                for name in type_sets:
                    if f"{name} types" in line:
                        tmp_content[i] = f"     {len(type_sets[name])}  {name} types\n"
                        break
                else:
                    if "Atoms" in line:
                        break

        output.writelines(tmp_content)


def generate_orca_input(input_file: str, xyz_filename: str, output_file: str) -> None:
    """
    Generate an ORCA input file using the specified input and .xyz file.

    Parameters
    ----------
    input_file : str
        Path to the input file.
    xyz_filename : str
        Path to the .xyz file.
    output_file : str
        Path to the output file.
    """
    # Check if the file exists and has the correct extension
    if not os.path.isfile(xyz_filename):
        print("Error: Please provide a valid .xyz file.")
        return
    if not os.path.isfile(input_file):
        print("Error: Please provide a valid input file.")
        return

    if not xyz_filename.endswith(".xyz"):
        xyz_filename = xyz_filename.replace(".xyz", "")
        xyz_filename += ".xyz"

    with open(input_file, "r") as file:
        lines = file.readlines()
    # Modify the line containing 'xyz_filename'
    success = False
    for i, line in enumerate(lines):
        if "*xyzfile" in line:
            success = True
            parts = line.split()
            # Replace the last part with the new filename
            if len(parts) == 4:
                parts[-1] = xyz_filename
                lines[i] = " ".join(parts) + "\n"
            elif len(parts) == 3:
                parts.append(f" {xyz_filename}")
            break
    if not success:
        print("Error: Please provide a valid input file containing *xyzfile signature")
        return
    # Write the modified lines to a new output file
    with open(output_file, "w") as file:
        file.writelines(lines)


def modify_lammps_data(data_file: str, xyz_file: str) -> None:
    """
    Modify the LAMMPS data file based on coordinates from the .xyz file.

    Parameters
    ----------
    data_file : str
        Path to the LAMMPS data file.
    xyz_file : str
        Path to the .xyz file.
    """
    # Read the xyz file
    with open(xyz_file, "r") as f:
        xyz_lines = f.readlines()

    # Extract coordinates from the xyz file
    num_atoms = int(xyz_lines[0].strip())
    coordinates = []
    for line in xyz_lines[2 : 2 + num_atoms]:
        parts = line.split()
        coordinates.append(parts[1:])  # Append x, y, z coordinates

    # Read and modify the .data file
    with open(data_file, "r") as f:
        data_lines = f.readlines()

    in_atoms_section = False
    atom_count = 0
    for i, line in enumerate(data_lines):
        if "Atoms" in line:
            in_atoms_section = True
        if in_atoms_section:
            # Modify lines in the "Atoms" section
            parts = line.split()
            # Replace the last three elements with new coordinates
            parts[-3:] = coordinates[atom_count]
            data_lines[i] = " ".join(parts) + "\n"
            atom_count += 1
            if atom_count >= num_atoms:
                break

    output_name = xyz_file.replace(".xyz", ".data")
    # Write the modified content to a new output file
    with open(output_name, "w") as f:
        f.writelines(data_lines)


def cut_ff(
    input_filename: str = "system.in.settings",
    style: str = "lj/cut/coul/long",
    units: str = "metal",
) -> None:
    """
    Cut force field data from the input file and save it into separate coefficient files.

    Parameters
    ----------
    input_filename : str
        Path to the input force field file.
    style : str
        Style for pair coefficients.
    units : str
        Units used in the force field file.
    """
    if units == "real":
        scaling_factor = 1
    elif units == "metal":
        scaling_factor = 0.0433641  # 1 kcal/mol = 0.0433641 eV
    # Get the directory of the input file
    input_dir = os.path.dirname(os.path.abspath(input_filename))

    # Define the output file paths in the same directory as the input file
    pair_file_path = os.path.join(input_dir, "pair_coefs.ff")
    other_file_path = os.path.join(input_dir, "other_coefs.ff")

    # Open the input file and two output files
    with (
        open(input_filename, "r") as ff_file,
        open(pair_file_path, "w") as pair_file,
        open(other_file_path, "w") as other_file,
        open(f"{input_filename}_{units}", "w") as ff_file_units,
    ):
        # Process the input file
        for line in ff_file:
            if "pair_coeff" in line:
                parsed = line.split(" ")
                start, end = parsed[:3], parsed[3:]
                end[0] = str(float(end[0]) * scaling_factor)
                pair_file.write(f"{' '.join(start)} {style} {' '.join(end)}")
                ff_file_units.write(f"{' '.join(start)} {' '.join(end)}")
            elif (
                ("bond_coeff" in line)
                or ("angle_coeff" in line)
                or ("improper_coeff" in line)
            ):
                parsed = line.split(" ")
                parsed[2] = str(float(parsed[2]) * scaling_factor)
                other_file.write(f"{' '.join(parsed)}")
                ff_file_units.write(f"{' '.join(parsed)}")
            elif "dihedral_coeff" in line:
                parsed = line.split(" ")
                parsed[2:6] = [str(x * scaling_factor) for x in map(float, parsed[2:6])]
                other_file.write(f"{' '.join(parsed)}\n")
                ff_file_units.write(f"{' '.join(parsed)}\n")
            else:
                other_file.write(line)
                ff_file_units.write(line)


def cut_nbh(
    input_name: str,
    out_base: str,
    atom_coords: list[float],
    radius: float,
    cut: bool,
    save_xyz: bool = False,
    save_pdb: bool = False,
) -> None:
    """
    Cut neighborhood atoms from a molecular structure file.

    Parameters
    ----------
    input_name : str
        Path to the input file.
    out_base : str
        Base name for the output files.
    atom_coords : list of float
        List of three floats representing (x, y, z) coordinates.
    radius : float
        Radius for selecting atoms.
    cut : bool
        Whether to cut the molecule.
    save_xyz : bool
        Whether to save .xyz output.
    save_pdb : bool
        Whether to save .pdb output.
    """

    if len(atom_coords) != 3:
        raise ValueError(
            "atom_coords must be a list of exactly three float values [x, y, z]."
        )

    out_base = os.path.splitext(out_base)[0]  # Remove extension if any
    input_ext = os.path.splitext(input_name)[1].lower()  # Get file extension

    # Initialize molecule object
    mol = molecule.Molecule()
    mol.input_geometry = input_name

    # Determine the correct backend based on input file type
    backend_mapping = {".xyz": "cclib", ".dump": "MDAnalysis", ".cfg": "internal"}

    backend = backend_mapping.get(input_ext)
    if backend is None and input_ext not in backend_mapping:
        raise ValueError(f"Unsupported input file format: {input_ext}")

    mol.build_geometry(backend=backend)

    # Process molecule
    mol.build_graph()
    mol.build_structure(cut_molecule=cut)

    # Select neighborhood atoms
    new_mol, idxs = mol.select(mol.select_r(atom_coords, radius))

    # Save outputs based on input type
    save_methods = {
        ".xyz": lambda: new_mol.save_xyz(out_base + ".xyz"),
        ".dump": lambda: mol.save_dump(out_base + ".dump", mol.input_geometry, idxs),
        ".cfg": lambda: mol.save_cfg(out_base + ".cfg", mol.input_geometry, idxs),
    }

    if input_ext in save_methods:
        save_methods[input_ext]()

    if save_xyz and not input_ext == ".xyz":
        new_mol.save_xyz(out_base + ".xyz")
    if save_pdb:
        new_mol.save_pdb(out_base + ".pdb")

    print(f"Processed {input_name} and saved results to {out_base}.")


def cut_random_nbhs(
    input_name: str,
    out_base: str,
    atom_coords: list[float],
    radius: float,
    cut: bool,
    save_xyz: bool = False,
    save_pdb: bool = False,
) -> None:
    pass


def cut_dump(
    input_file: str, output_prefix: str = "frame_", starting_id: int = 0
) -> int:
    """
    Splits a LAMMPS dump file into separate files, each containing a single timestep.

    Parameters
    ----------
    input_file : str
        Path to the input dump file.
    output_prefix : str, default="frame\_"
        Prefix for output files.

    Returns
    -------
    int
        Last frame id.
    """
    frame = starting_id
    outfile = None

    with open(input_file, "r", encoding="utf-8") as dump:
        lines = dump.readlines()

    num_atoms = int(lines[3])

    for line in lines:
        if "ITEM: TIMESTEP" in line:
            # Close the previous file if it was open
            if outfile:
                outfile.close()

            # Create a new file for the current frame
            output_filename = f"{output_prefix}{frame:04d}.dump"
            outfile = open(output_filename, "w", encoding="utf-8")
            frame += 1

        # Write the current line to the corresponding file
        if outfile:
            outfile.write(line)

    # Close the last open file
    if outfile:
        outfile.close()

    print(f"Splitting completed: {frame} files created.")

    return frame - 1


def filter_by_grade(input_file: str, output_file: str, threshold: float = 2.0) -> None:
    """
    Filter configurations based on a minimum grade threshold and save them to the output file.

    Parameters
    ----------
    input_file : str
        Path to the input file.
    output_file : str
        Path to the output file where filtered results will be stored.
    n : float, default: 2.0
        Minimum grade threshold; only configurations with a grade >= n will be included.
    """
    with open(input_file, "r", encoding="utf-8") as infile, open(
        output_file, "w", encoding="utf-8", newline="\n"
    ) as outfile:
        inside_cfg = False
        current_cfg = []
        current_grade = None

        for line in infile:
            line = line.strip()

            if line == "BEGIN_CFG":
                inside_cfg = True
                current_cfg = [line]
                current_grade = None
            elif line == "END_CFG":
                if (
                    inside_cfg
                    and current_grade is not None
                    and current_grade >= threshold
                ):
                    current_cfg.append(line)
                    outfile.write("\n".join(current_cfg) + "\n\n")
                inside_cfg = False
                current_cfg = []  # Reset buffer
            elif inside_cfg:
                if line.startswith("Feature   MV_grade"):
                    current_grade = float(line.split()[-1])
                    if current_grade < threshold:
                        inside_cfg = False  # Discard cfg
                        current_cfg = []  # Clear buffer
                        continue
                current_cfg.append(line)

    print(
        f"Configurations with grade greater than {threshold} have been extracted to {output_file}."
    )


def check_cfg_ratio(
    new_train_file: str, current_train_file: str, fraction: float = 0.5
) -> bool:
    """
    Compare the number of 'BEGIN_CFG' occurrences in two .cfg files and return True if
    the count in new_train is at least fraction times the count in current_train.

    Parameters
    ----------
    new_train_file : str
        Path to the new training .cfg file.
    current_train_file : str
        Path to the current training .cfg file.
    fraction : float, optional, default: 0.5
        The fraction threshold.

    Returns
    -------
    bool
        True if the count of 'BEGIN_CFG' in new_train is >= fraction * count in current_train, otherwise False.
    """

    def count_begin_cfg(file_path: str) -> int:
        with open(file_path, "r", encoding="utf-8") as f:
            return sum(1 for line in f if line.strip() == "BEGIN_CFG")

    new_count = count_begin_cfg(new_train_file)
    current_count = count_begin_cfg(current_train_file)

    return new_count >= fraction * current_count
