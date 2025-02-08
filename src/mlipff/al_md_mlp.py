import os
import time
import subprocess
from pathlib import Path


# Function to run a shell command
# 'shell=True' allows commands as strings, similar to writing them directly in a shell
# 'check=True' raises an error if the command fails
def run_command(command):
    subprocess.run(command, shell=True, check=True)


# Function to check if a file exists using Pathlib
# Path(filepath).exists() returns True if the file exists, False otherwise
def is_file_exists(filepath):
    return Path(filepath).exists()


# Function to batch remove files using glob patterns
# Path().glob(file) allows for pattern matching, similar to using '*' in the shell
def batch_remove(files):
    for file in files:
        for path in Path().glob(file):
            path.unlink(
                missing_ok=True
            )  # unlink removes the file, missing_ok avoids errors if the file doesn't exist


n_cores = 16  # Number of cores used, adjust accordingly

# Cleanup files before starting to ensure no leftover data
batch_remove(
    ["preselected.cfg.*", "preselected.cfg", "selected.cfg", "nbh.cfg", "*.txt"]
)

# Main loop to continuously run the workflow
while True:
    # Backup the current almtp file and create necessary files
    run_command("cp curr.almtp pot_save/")
    Path("preselected.cfg").touch()  # Creates an empty file named 'preselected.cfg'

    # Run LAMMPS using sbatch
    is_lammps_finished = False
    run_command(
        f"sbatch -J lammps -o lammps.out -e lammps.err -N 1 -n {n_cores} ./run_lammps.sh"
    )

    # Poll until LAMMPS is finished
    while not is_lammps_finished:
        time.sleep(5)  # Sleep for 5 seconds between checks
        is_lammps_finished = is_file_exists("is_lammps_finished.txt")
        print(int(is_lammps_finished))  # Print status (0 or 1)

    time.sleep(5)  # Extra sleep to ensure processes have completed

    # Concatenate all preselected configuration files into one
    run_command("cat preselected.cfg.* >> preselected.cfg")
    batch_remove(["preselected.cfg.*"])  # Remove intermediate files

    # Count the number of preselected configurations
    n_preselected = sum(1 for line in Path("preselected.cfg").open() if "BEGIN" in line)
    if n_preselected > 0:
        print("selection")

        # Run selection job using sbatch
        is_selection_finished = False
        run_command(
            "sbatch -J selection -o selection.out -e selection.err -N 1 -n 8 ./select.sh"
        )

        # Poll until selection job is finished
        while not is_selection_finished:
            time.sleep(5)
            is_selection_finished = is_file_exists("is_selection_finished.txt")
            print(int(is_selection_finished))

        time.sleep(1)

        # Check if selected.cfg exists; if not, exit the loop
        if not is_file_exists("selected.cfg"):
            print("selected.cfg does not exist")
            break

        # Process selected configurations
        batch_remove(["preselected.cfg"])
        run_command("./mlipff cut_nbh --input=selected.cfg --output=nbh.cfg --cutoff=5")
        batch_remove(["selected.cfg"])

        # Count the number of POSCAR configurations
        n_orca = sum(1 for line in Path("nbh.cfg").open() if "BEGIN" in line)
        n_orca_prev = n_orca - 1

        # Organize Orca directories and files
        for i in range(n_orca):
            orca_dir = Path(f"ORCA/{i}")
            orca_dir.mkdir(
                parents=True, exist_ok=True
            )  # Create directory if it doesn't exist
            run_command(f"cp ORCA_input_{i}.inp ORCA/{i}/")
            run_command(f"cp ORCA_input_{i}.xyz ORCA/{i}/")
            # Copy necessary files into each Orca directory

        # Submit QM jobs for each configuration
        for i in range(n_orca):
            os.chdir(f"ORCA/{i}")
            run_command("sbatch ./run_orca.sh")
            os.chdir(
                Path.cwd().parent.parent
            )  # Navigate back to the original working directory

        # Wait for QM jobs to finish with adaptive sleep
        is_orca_finished = False
        sleep_interval = 5
        while not is_orca_finished:
            time.sleep(sleep_interval)
            # Count how many Orca jobs have finished by checking for 'is_orca_finished.txt'
            n_qm_finished = len(list(Path("ORCA").glob("*/is_orca_finished.txt")))
            is_orca_finished = n_orca == n_qm_finished
            print(int(is_orca_finished))
            # Increase sleep interval up to a maximum of 60 seconds
            sleep_interval = min(sleep_interval * 2, 60)

        time.sleep(1)

        # Process Orca log files to generate calculated.cfg
        for i in range(n_orca):
            run_command(
                f"./mlipff convert --input=ORCA/{i}/out.log --output=ORCA/{i}/calculated.cfg --from=orcadump --to=cfg"
            )
            with Path(f"ORCA/{i}/calculated.cfg").open("r") as cfg_file:
                with Path("train.cfg").open("a") as train_file:
                    train_file.write(
                        cfg_file.read()
                    )  # Append contents of each calculated.cfg to train.cfg
            # batch_remove([f"ORCA/{i}"])  # Remove the ORCA directory after processing

        # Run training job using sbatch
        is_training_finished = False
        run_command(
            "sbatch -J train -o train.out -e train.err -N 1 --exclusive  ./train.sh --time=3:00:00"
        )

        # Poll until training is finished
        while not is_training_finished:
            time.sleep(5)
            is_training_finished = is_file_exists("is_training_finished.txt")
            print(int(is_training_finished))

        # Cleanup unnecessary files after training
        time.sleep(5)
        batch_remove(["nbh.cfg", "*.txt"])
    else:
        break
