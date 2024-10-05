class Configuration:
    """
    A class to represent a single configuration of atoms.

    Attributes
    ----------
    size : int
        The number of atoms in the configuration.
    atom_data : list of dict
        The list of atom data with properties id, type, x, y, z, fx, fy, fz.
    supercell : list of float, optional
        The supercell matrix describing the configuration.
    energy : float, optional
        The energy of the configuration.
    plus_stress : dict, optional
        Stress information of the configuration.
    features : dict, optional
        Additional features of the configuration.
    """

    def __init__(
        self,
        size=None,
        atom_data=None,
        supercell=None,
        energy=None,
        plus_stress=None,
        features=None,
    ):
        self.size = size
        self.atom_data = atom_data
        self.supercell = supercell or []
        self.energy = energy or 0
        self.plus_stress = plus_stress or {}
        self.features = features or {}

    def read_cfg(self, file_path):
        """Reads a single configuration from a .cfg file."""
        with open(file_path, "r") as file:
            atom_data = []
            reading_atoms = False
            file_iter = iter(file)

            for line in file_iter:
                line = line.strip()

                if "BEGIN_CFG" in line:
                    continue

                if "Size" in line:
                    self.size = int(next(file).strip())
                    continue
                elif "Supercell" in line:
                    self.supercell = []
                    for _ in range(3):
                        line = next(file).strip()
                        if "AtomData:" in line:
                            reading_atoms = True
                            atom_headers = line.split()[1:]
                            break
                        else:
                            self.supercell.extend(map(float, line.split()))
                    continue
                elif "AtomData:" in line:
                    reading_atoms = True
                    atom_headers = line.split()[1:]
                    continue
                elif reading_atoms:
                    if line == "":
                        reading_atoms = False
                    else:
                        atom_values = line.split()
                        atom = {}
                        for i, value in enumerate(atom_values):
                            if atom_headers[i] in ["id", "type"]:
                                atom[atom_headers[i]] = int(value)
                            else:
                                atom[atom_headers[i]] = float(value)
                        atom_data.append(atom)
                    if len(atom_data) == self.size:
                        reading_atoms = False
                        self.atom_data = atom_data
                        continue
                elif "Energy" in line:
                    self.energy = float(next(file).strip())
                elif "PlusStress" in line:
                    stress_headers = line.split()[1:]
                    stress_values = list(map(float, next(file).strip().split()))
                    self.plus_stress = dict(zip(stress_headers, stress_values))
                elif "Feature" in line:
                    feature_parts = line.split(maxsplit=2)
                    if len(feature_parts) == 3:
                        feature_name = feature_parts[1]
                        feature_value = feature_parts[2]
                        self.features[feature_name] = feature_value

    def write_cfg(self, file_path):
        """Writes the configuration to a .cfg file."""
        with open(file_path, "w") as file:
            file.write("BEGIN_CFG\n")
            file.write(f"Size\n{self.size}\n")
            if self.supercell:
                file.write("Supercell\n")
                for i in range(0, len(self.supercell), 3):
                    file.write(f"{' '.join(map(str, self.supercell[i:i+3]))}\n")
            if self.atom_data:
                headers = " ".join(self.atom_data[0].keys())
                file.write(f"AtomData: {headers}\n")
                for atom in self.atom_data:
                    file.write(" ".join(map(str, atom.values())) + "\n")
            if self.energy is not None:
                file.write(f"Energy\n{self.energy}\n")
            if self.plus_stress:
                file.write("PlusStress: " + " ".join(self.plus_stress.keys()) + "\n")
                file.write(" ".join(map(str, self.plus_stress.values())) + "\n")
            if self.features:
                for name, value in self.features.items():
                    file.write(f"Feature {name} {value}\n")
            file.write("END_CFG\n\n")

    def read_dump(self, file_path, include_efs, replace_types):
        """
        Reads a LAMMPS .dump file and populates the configuration attributes.

        Parameters
        ----------
        file_path : str
            The path to the .dump file to read from.
        include_efs : bool
            Whether to include energy, forces, and stress data.
        """
        type_replacement_dict = {}
        if replace_types:
            with open(replace_types, "r") as file:
                for line in file:
                    old_type, new_type = line.strip().split()
                    type_replacement_dict[int(old_type)] = int(new_type)

        with open(file_path, "r") as lmp_dump:
            lines = lmp_dump.readlines()

        lattice = [0.0, 0.0, 0.0]
        data = []
        in_atoms_section = False
        in_box_section = False

        for line in lines:
            if "ITEM: NUMBER" in line:
                self.size = int(lines[lines.index(line) + 1].strip())
            elif "ITEM: BOX" in line:
                in_box_section = True
            elif in_box_section:
                lattice[lines.index(line) % 3] = float(line.split()[1])
                if lines.index(line) % 3 == 2:
                    in_box_section = False
                    self.supercell = [
                        lattice[0],
                        0,
                        0,
                        0,
                        lattice[1],
                        0,
                        0,
                        0,
                        lattice[2],
                    ]
            elif "ITEM: ATOMS" in line:
                in_atoms_section = True
                coordinates_start_index = line.split().index("x") - 2
                if "fx" in line:
                    include_efs = True
                    forces_start_index = line.split().index("fx") - 2
                else:
                    include_efs = False
            elif in_atoms_section:
                atoms_data = list(map(float, line.split()))
                atom_type = int(atoms_data[1] - 1)
                if atom_type in type_replacement_dict:
                    atom_type = int(type_replacement_dict[atom_type + 1] - 1)
                pos_x, pos_y, pos_z = atoms_data[
                    coordinates_start_index : coordinates_start_index + 3
                ]
                f_x, f_y, f_z = (
                    atoms_data[forces_start_index : forces_start_index + 3]
                    if include_efs
                    else (0, 0, 0)
                )
                self.energy += (
                    atoms_data[forces_start_index + 3] if include_efs else None
                )
                atom = {
                    "id": len(data) + 1,
                    "type": atom_type,
                    "cartes_x": pos_x,
                    "cartes_y": pos_y,
                    "cartes_z": pos_z,
                    "fx": f_x,
                    "fy": f_y,
                    "fz": f_z,
                }
                data.append(atom)

        self.atom_data = data

    def read_orca_dump(self, filename, replace_types):
        """
        Reads ORCA .log and LAMMPS .dump files, creating substracted configuration
        QM - MM

        Parameters
        ----------
        filename : str
            The base filename (without extension) of the .log and .dump files.
        """

        if len(filename.split(".")) != 1:
            filename = filename.split(".")[0]
        logname = f"{filename}.log"
        dumpname = f"{filename}.dump"

        # Constants
        Eh = 27.2113834
        Bohr = 0.5291772083

        # Read initial configuration from .dump file
        mm_config = Configuration()
        mm_config.read_dump(dumpname, True, replace_types)

        # Create a new configuration with the same atomic positions but different forces/energy
        qm_config = Configuration(
            size=mm_config.size,
            atom_data=[
                {
                    "id": atom["id"],
                    "type": atom["type"],
                    "cartes_x": atom["cartes_x"],
                    "cartes_y": atom["cartes_y"],
                    "cartes_z": atom["cartes_z"],
                    "fx": atom["fx"],
                    "fy": atom["fy"],
                    "fz": atom["fz"],
                }
                for atom in mm_config.atom_data
            ],
            supercell=mm_config.supercell,
            energy=mm_config.energy,
        )

        # Read forces and energy from .log file
        with open(logname, "r") as logfile:
            lines = logfile.readlines()

        try:
            energy_line = next(
                line for line in lines if "FINAL SINGLE POINT ENERGY" in line
            )
            qm_energy = float(energy_line.split()[4]) * Eh - mm_config.energy
        except StopIteration:
            raise ValueError("Energy not found in ORCA log file")

        qm_config.energy = qm_energy

        try:
            grad_start = lines.index("CARTESIAN GRADIENT\n") + 3
            for atom, line in zip(qm_config.atom_data, lines[grad_start:]):
                force_data = list(map(float, line.split()[3:6]))
                atom["fx"] = -force_data[0] * Eh / Bohr - atom["fx"]
                atom["fy"] = -force_data[1] * Eh / Bohr - atom["fy"]
                atom["fz"] = -force_data[2] * Eh / Bohr - atom["fz"]
        except (ValueError, IndexError):
            raise ValueError("Gradient data not found or incomplete in ORCA log file")

        # Replace the original configuration's atom data and energy with updated values
        self.size = qm_config.size
        self.atom_data = qm_config.atom_data
        self.supercell = []
        self.energy = qm_config.energy

    def __sub__(self, other):
        """Subtracts the fields fx, fy, fz, Energy, and PlusStress from another Configuration object."""
        new_atom_data = []
        for atom1, atom2 in zip(self.atom_data, other.atom_data):
            new_atom = {}
            for key in atom1:
                if key in ["fx", "fy", "fz"]:
                    new_atom[key] = atom1.get(key, 0.0) - atom2.get(key, 0.0)
                else:
                    new_atom[key] = atom1[key]
            new_atom_data.append(new_atom)

        new_plus_stress = {
            key: self.plus_stress.get(key, 0.0) - other.plus_stress.get(key, 0.0)
            for key in self.plus_stress
        }

        return Configuration(
            size=self.size,
            atom_data=new_atom_data,
            supercell=self.supercell,
            energy=(
                (self.energy - other.energy)
                if self.energy is not None and other.energy is not None
                else None
            ),
            plus_stress=new_plus_stress,
            features=self.features,
        )
