import random
import os
import shutil
from mlipff_configuration import Configuration
from ocellar import molecule


def convert(from_format, to_format, input_file, output_file):
    config = Configuration()

    if from_format == "dump":
        config.read_dump(input_file)
    elif from_format == "orca":
        config.read_orca_dump(input_file)
    else:
        print(f"Unsupported input format: {from_format}")
        return

    if to_format == "cfg":
        config.write_cfg(output_file)
    else:
        print(f"Unsupported output format: {to_format}")


def cut_data(data_file, dump_file, output_file):
    tmp_num = random.randint(100000, 999999)
    tmp_dir = f"data_tmp_{tmp_num}"
    os.mkdir(tmp_dir)
    tmp1_path = os.path.join(tmp_dir, "tmp1.txt")
    tmp2_path = os.path.join(tmp_dir, "tmp2.txt")
    tmp3_path = os.path.join(tmp_dir, "tmp3.txt")

    with (
        open(data_file, "r") as datafile,
        open(dump_file, "r") as dumpfile,
        open(tmp1_path, "w") as tmp1,
        open(tmp2_path, "w") as tmp2,
        open(tmp3_path, "w") as tmp3,
        open(output_file, "w") as output,
    ):
        ids = set()
        read_number = False
        read_atoms = False

        for line in dumpfile:
            if "ITEM: NUMBER" in line:
                read_number = True
                continue
            if read_number:
                atom_number = int(line.strip())
                read_number = False
            if "ITEM: ATOMS" in line:
                read_atoms = True
                continue
            if read_atoms:
                ids.add(int(line.split()[0]))
                if len(ids) == atom_number:
                    break

        parsarg = ""
        line_counter = 0
        atoms = bonds = angles = dihedrals = impropers = 0

        for line in datafile:
            if line_counter <= 23:
                tmp1.write(line)
                line_counter += 1
                continue
            if line == "\n":
                tmp2.write(line)
                continue
            if "Atoms" in line:
                parsarg = "atoms"
                tmp2.write(line)
                continue
            if "Bonds" in line:
                parsarg = "bonds"
                tmp2.write(line)
                continue
            if "Angles" in line:
                parsarg = "angles"
                tmp2.write(line)
                continue
            if "Dihedrals" in line:
                parsarg = "dihedrals"
                tmp2.write(line)
                continue
            if "Impropers" in line:
                parsarg = "impropers"
                tmp2.write(line)
                continue
            if parsarg == "atoms":
                if int(line.split()[0]) in ids:
                    atoms += 1
                    tmp2.write(line)
            if parsarg == "bonds":
                if int(line.split()[2]) in ids or int(line.split()[3]) in ids:
                    bonds += 1
                    tmp2.write(line)
            if parsarg == "angles":
                if (
                    int(line.split()[2]) in ids
                    or int(line.split()[3]) in ids
                    or int(line.split()[4]) in ids
                ):
                    angles += 1
                    tmp2.write(line)
            if parsarg == "dihedrals":
                if (
                    int(line.split()[2]) in ids
                    or int(line.split()[3]) in ids
                    or int(line.split()[4]) in ids
                    or int(line.split()[5]) in ids
                ):
                    dihedrals += 1
                    tmp2.write(line)
            if parsarg == "impropers":
                if (
                    int(line.split()[2]) in ids
                    or int(line.split()[3]) in ids
                    or int(line.split()[4]) in ids
                    or int(line.split()[5]) in ids
                ):
                    impropers += 1
                    tmp2.write(line)

        with open(tmp1_path, "r") as tmp1:
            for line in tmp1:
                if "atoms" in line:
                    tmp3.write(f"     {atoms}  atoms\n")
                elif "bonds" in line:
                    tmp3.write(f"     {bonds}  bonds\n")
                elif "angles" in line:
                    tmp3.write(f"     {angles}  angles\n")
                elif "dihedrals" in line:
                    tmp3.write(f"     {dihedrals}  dihedrals\n")
                elif "impropers" in line:
                    tmp3.write(f"     {impropers}  impropers\n")
                else:
                    tmp3.write(line)

        shutil.copyfileobj(open(tmp3_path, "r"), output)
        shutil.copyfileobj(open(tmp2_path, "r"), output)

    shutil.rmtree(tmp_dir)


def generate_orca_input(input_file, xyz_filename: str):
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
    out_filename = xyz_filename.replace(".xyz", ".inp")
    # Write the modified lines to a new output file
    with open(out_filename, "w") as file:
        file.writelines(lines)


def modify_lammps_data(data_file, xyz_file):
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
            continue
        if in_atoms_section and line.strip() == "":
            # End of "Atoms" section
            in_atoms_section = False
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


def cut_ff(input_filename="system.in.settings", style="lj/cut/coul/long", units = 'metal'):
    if units == 'real':
        scaling_factor = 1
    elif units == 'metal':
        scaling_factor = 0.0433641 # 1 kcal/mol = 0.0433641 eV
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
    ):

        # Process the input file
        for line in ff_file:
            if "pair_coeff" in line:
                parsed = line.split(" ")
                start, end = parsed[:3], parsed[3:]
                end[0] *= scaling_factor
                pair_file.write(f"{' '.join(start)} {style} {' '.join(end)}")
            elif ("bond_coeff" in line) or ("angle_coeff" in line) or ("improper_coeff" in line):
                parsed = line.split(" ")
                parsed[2] *= scaling_factor
                pair_file.write(f"{' '.join(parsed)}")
            elif ("dihedral_coeff" in line):
                parsed = line.split(" ")
                parsed[2:6] *= scaling_factor
                pair_file.write(f"{' '.join(parsed)}")
            else:
                other_file.write(line)


def cut_nbh(input_name, out_base, atom_coords, radius, cut):
    out_dump = str(out_base) + ".dump"
    out_xyz = str(out_base) + ".xyz"
    # from dump
    mol = molecule.Molecule()
    mol.input_geometry = input_name
    mol.build_geometry(
        backend="MDAnalysis"
    )  # Build geometry from input file with MDAnalysis backend
    mol.build_graph()
    mol.build_structure(cut_molecule=cut)

    new_mol, idxs = mol.select(
        mol.select_r(atom_coords, radius)
    )  # Select atom idxs with sphere center and radius, then build a new Molecule and idxs

    mol.save_dump(
        out_dump, mol.input_geometry, idxs
    )  # Save dump file using original dump and idxs

    new_mol.build_graph()
    new_mol.build_structure(cut_molecule=cut)
    new_mol.save_xyz(out_xyz)
