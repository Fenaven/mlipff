# mlipff
MLIP-FF Project
Description
MLIP-FF is a Python-based tool for using MTP along with classical force-fields inside LAMMPS with additional usage of MLIP. It supports various operations on molecular simulation data, including file conversion, partitioning force fields, and generating input files for quantum chemistry software such as Orca.

The project is set up using PDM (Python Development Master) for dependency management and can be easily extended or adapted for other use cases.

Features
Convert Files: Convert files between different formats (e.g., from Orca or dump to cfg).
Cut Data: Partition and manipulate system data using dump files.
Generate Orca Input: Automatically generate input files for Orca from XYZ coordinate files.
Force Field Manipulation: Partition force field files based on styles.
Dependencies
ocellar: A package for molecular structure manipulation, installed from GitHub.
Other Python libraries (listed in the pyproject.toml file).
Installation
To install the MLIPFF project and its dependencies, you can use PDM.

Clone the repository:

bash
Copy code
git clone https://github.com/Fenaven/mlipff.git
cd mlipff
Install the dependencies using PDM:

bash
Copy code
pdm install
Usage
The main script for interacting with the tool is mlipff.py. You can run it with different commands to perform various operations.

Example Commands:
Convert files:

bash
Copy code
python mlipff.py convert --from dump --to cfg --input input_file.dump --output output_file.cfg
Cut data:

bash
Copy code
python mlipff.py cut_data --data system.data --nbh nbh.dump --output output.data
Generate Orca input:

bash
Copy code
python mlipff.py create_input --xyz structure.xyz
Cut force field file:

bash
Copy code
python mlipff.py cut_ff --input system.in.settings --style lj/cut/coul/cut
Project Structure
mlipff.py: The main entry point of the project.
mlipff_configuration.py: Handles the configuration of molecular systems.
utils.py: Provides utility functions for file conversion, data manipulation, and force field partitioning.
tests/: (Optional) Placeholder for unit tests.