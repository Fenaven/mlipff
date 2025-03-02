import os
from typing import Optional, List, Dict, Union


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
        The simulation box bounds stored as [xmin, xmax, ymin, ymax, zmin, zmax].
    energy : float, optional
        The energy of the configuration.
    plus_stress : dict, optional
        Stress information of the configuration.
    features : dict, optional
        Additional features of the configuration.
    """

    def __init__(
        self,
        size: Optional[int] = None,
        atom_data: Optional[List[Dict[str, Union[int, float]]]] = None,
        supercell: Optional[List[float]] = None,
        energy: Optional[float] = None,
        plus_stress: Optional[Dict[str, float]] = None,
        features: Optional[Dict[str, str]] = None,
    ) -> None:
        """
        Initializes the Configuration object with optional attributes.

        Parameters
        ----------
        size : Optional[int]
            The number of atoms in the configuration.
        atom_data : Optional[List[Dict[str, Union[int, float]]]]
            The list of atom data.
        supercell : Optional[List[float]]
            The supercell bounds stored as [xmin, xmax, ymin, ymax, zmin, zmax].
        energy : Optional[float]
            The energy of the configuration.
        plus_stress : Optional[Dict[str, float]]
            Stress information of the configuration.
        features : Optional[Dict[str, str]]
            Additional features of the configuration.
        """
        self.size = size
        self.atom_data = atom_data
        self.supercell = supercell if supercell is not None else []
        self.energy = energy or 0
        self.plus_stress = plus_stress or {}
        self.features = features or {}

    def read_cfg(self, file_path: str, store_grade: Optional[int]) -> None:
        """
        Reads a single configuration from a .cfg file.

        Parameters
        ----------
        file_path : str
            The path to the .cfg file to read from.
        """
        self.atom_data = []
        self.supercell = []
        self.energy = None
        self.plus_stress = {}
        self.features = {}

        with open(file_path, "r") as file:
            lines = [
                line.strip() for line in file
            ]  # Read all lines and strip whitespace

        for i, line in enumerate(lines):
            if line == "Size":
                self.size = int(lines[i + 1])

            elif line == "Supercell":
                self.supercell = [0, 0, 0, 0, 0, 0]
                for axis in range(3):
                    values = list(map(float, lines[axis + 1]))
                    self.supercell[axis * 2 + 1] = values[
                        axis
                    ]  # Store xmax, ymax, zmax

            elif line.startswith("AtomData:"):
                atom_headers = line.split()[1:]
                self.atom_data = []
                for j in range(self.size):
                    values = lines[i + 1 + j].split()  # Read next self.size lines
                    atom = {
                        key: int(val) if key in ["id", "type"] else float(val)
                        for key, val in zip(atom_headers, values)
                    }
                    self.atom_data.append(atom)

                i += self.size

            elif line == "Energy":
                self.energy = float(lines[i + 1])

            elif line.startswith("PlusStress"):
                stress_headers = line.split()[1:]
                stress_values = map(float, lines[i + 1].split())
                self.plus_stress = dict(zip(stress_headers, stress_values))

            elif line.startswith("Feature"):
                feature_parts = line.split(maxsplit=2)
                if len(feature_parts) == 3:
                    self.features[feature_parts[1]] = feature_parts[2]

    def old_read_cfg(self, file_path: str) -> None:
        """
        Reads a single configuration from a .cfg file.

        Parameters
        ----------
        file_path : str
            The path to the .cfg file to read from.
        """
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
                    self.supercell = [0, 0, 0, 0, 0, 0]  # Initialize
                    for i, axis in enumerate(["x", "y", "z"]):
                        line = next(file).strip()
                        values = list(map(float, line.split()))
                        self.supercell[i * 2 + 1] = values[i]  # Store xmax, ymax, zmax
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

    def write_cfg(self, file_path: str) -> None:
        """
        Writes the configuration to a .cfg file.

        Parameters
        ----------
        file_path : str
            The path to the .cfg file to write to.
        """
        with open(file_path, "w") as file:
            file.write("BEGIN_CFG\n")
            file.write(f"Size\n{self.size}\n")
            if self.supercell:
                file.write("Supercell\n")
                file.write(f"{self.supercell[1] - self.supercell[0]} 0 0\n")
                file.write(f"0 {self.supercell[3] - self.supercell[2]} 0\n")
                file.write(f"0 0 {self.supercell[5] - self.supercell[4]}\n")
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

    def read_dump(
        self, file_path: str, include_efs: bool, replace_types: Optional[str]
    ) -> None:
        """
        Reads a LAMMPS .dump file and populates the configuration attributes.
            Columns of dump file should be:
            id type mass x y z fx fy fz c_pe
            id and type should be always on 1 and 2 place, position of xyz and forces+energy may vary,
            but the order x y z should be followed, as well as fx fy fz c_pe

        Parameters
        ----------
        file_path : str
            The path to the .dump file to read from.
        include_efs : bool
            Whether to include energy, forces, and stress data.
        replace_types : Optional[str]
            The path to the file containing type replacement information.
        """
        type_replacement_dict = {}
        if replace_types:
            with open(replace_types, "r") as file:
                for line in file:
                    old_type, new_type = line.strip().split()
                    type_replacement_dict[int(old_type)] = int(new_type)

        with open(file_path, "r") as lmp_dump:
            lines = lmp_dump.readlines()

        data = []
        in_atoms_section = False

        for i, line in enumerate(lines):
            if "ITEM: NUMBER" in line:
                self.size = int(lines[i + 1].strip())

            elif "ITEM: BOX" in line:
                x_min, x_max = map(float, lines[i + 1].split())
                y_min, y_max = map(float, lines[i + 2].split())
                z_min, z_max = map(float, lines[i + 3].split())
                self.supercell = [x_min, x_max, y_min, y_max, z_min, z_max]

            elif "ITEM: ATOMS" in line:
                in_atoms_section = True
                headers = line.split()[2:]
                coordinates_start_index = headers.index("x")
                forces_start_index = headers.index("fx") if "fx" in headers else None

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

                self.energy += atoms_data[forces_start_index + 3] if include_efs else 0
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

    def read_orca(
        self, filename: str, substract: bool, replace_types: Optional[str]
    ) -> None:
        """
        Reads ORCA .log and LAMMPS .dump files, creating QM configuration.

        Parameters
        ----------
        filename : str
            The base filename (without extension) of the .log and .dump files.
        replace_types : Optional[str]
            The filename of a file which contains information for replacing atomic types.
        """

        # Extract the base name without extension
        base_filename = os.path.splitext(filename)[0]
        logname = f"{base_filename}.log"
        dumpname = f"{base_filename}.dump"

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
            qm_energy = float(energy_line.split()[4]) * Eh
            if substract:
                qm_energy -= mm_config.energy
        except StopIteration:
            raise ValueError("Energy not found in ORCA log file")

        qm_config.energy = qm_energy

        try:
            grad_start = lines.index("CARTESIAN GRADIENT\n") + 3
            for atom, line in zip(qm_config.atom_data, lines[grad_start:]):
                force_data = list(map(float, line.split()[3:6]))
                old_fxyz = [atom["fx"], atom["fy"], atom["fz"]]
                atom["fx"] = -force_data[0] * Eh / Bohr
                atom["fy"] = -force_data[1] * Eh / Bohr
                atom["fz"] = -force_data[2] * Eh / Bohr
                if substract:
                    atom["fx"] -= old_fxyz[0]
                    atom["fy"] -= old_fxyz[1]
                    atom["fz"] -= old_fxyz[2]
        except (ValueError, IndexError):
            raise ValueError("Gradient data not found or incomplete in ORCA log file")

        # Replace the original configuration's atom data and energy with updated values
        self.size = qm_config.size
        self.atom_data = qm_config.atom_data
        self.supercell = []
        self.energy = qm_config.energy

    def calculate_forces(self) -> float:
        forces_sum = 0
        for atom in self.atom_data:
            forces_sum += atom["fx"] ** 2 + atom["fy"] ** 2 + atom["fz"] ** 2
        return forces_sum

    def __sub__(self, other: "Configuration") -> "Configuration":
        """
        Subtracts the fields fx, fy, fz, Energy, and PlusStress from another Configuration object.

        Parameters
        ----------
        other : Configuration
            The other Configuration object to subtract from this one.

        Returns
        -------
        Configuration
            A new Configuration object with the subtracted fields.
        """
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
