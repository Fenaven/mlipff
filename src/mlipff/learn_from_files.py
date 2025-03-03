import os
import time
import random
import shutil
from pathlib import Path
from utils import (
    filter_by_grade,
    cut_data,
    convert,
    generate_orca_input,
    cut_nbh,
    check_cfg_ratio,
    cut_dump,
)
from mlipff_configuration import Configuration


def is_sphere_inside_cell(atom_coords, radius, box_bounds):
    """
    Checks if a sphere with a given center and radius is fully inside the simulation box.

    Parameters:
        atom_coords (tuple): (x, y, z) coordinates of the sphere center.
        radius (float): Radius of the sphere.
        box_bounds (list): List containing [x_min, x_max, y_min, y_max, z_min, z_max].

    Returns:
        bool: True if the sphere is inside the box, False otherwise.
    """
    x, y, z = atom_coords
    x_min, x_max, y_min, y_max, z_min, z_max = box_bounds

    return (
        x - radius >= x_min
        and x + radius <= x_max
        and y - radius >= y_min
        and y + radius <= y_max
        and z - radius >= z_min
        and z + radius <= z_max
    )


def preprocess_single_frame(dump_file, radius, cut):
    """
    Processes a single MD frame:
    - Reads the dump file.
    - Extracts box bounds and atom data.
    - Selects a random atom and checks if a sphere around it fits within the box.
    - If valid, performs a neighborhood cut.

    Parameters:
        dump_file (str): Path to the dump file.
        radius (float): Radius for the neighborhood selection.
        cut_molecules (bool): Whether to cut molecules in `cut_nbh`.

    Returns:
        bool: True if the cut was performed, False otherwise.
    """
    config = Configuration()

    try:
        config.read_dump(dump_file, include_efs=False, replace_types=None)
    except FileNotFoundError:
        print(f"Warning: {dump_file} not found.")
        return False

    if not config.atom_data:
        print(f"Skipping {dump_file}: No atom data available.")
        return False

    # Extract box bounds from the dump file
    box_bounds = config.supercell

    if box_bounds is None:
        print(f"Skipping {dump_file}: Box bounds not found.")
        return False

    # Select a random atom
    atom = random.choice(config.atom_data)
    atom_coords = (atom["cartes_x"], atom["cartes_y"], atom["cartes_z"])

    # Check if the sphere is fully inside the simulation box
    while not is_sphere_inside_cell(atom_coords, radius, box_bounds):
        atom = random.choice(config.atom_data)
        atom_coords = (atom["cartes_x"], atom["cartes_y"], atom["cartes_z"])

    output_name = "cut/" + dump_file.replace(".dump", "_cut")
    print(output_name)
    cut_nbh(
        input_name=dump_file,
        out_base=output_name,
        atom_coords=list(atom_coords),
        radius=radius,
        cut=cut,
        save_xyz=True,
    )
    print(f"Frame {dump_file}: Cut performed around atom {atom['id']} at {atom_coords}")
    return True


def preprocess_multiple_frames(num_frames, file_prefix, radius, cut_molecules):
    """
    Processes multiple MD frames:
    - Reads dump files based on user-defined prefix.
    - Calls `process_single_frame` for each file.
    - Logs the number of successful cuts.

    Parameters:
        num_frames (int): Number of frames to process.
        file_prefix (str): Prefix for the dump files.
        radius (float): Radius for the neighborhood selection.
        cut_molecules (bool): Whether to cut molecules in `cut_nbh`.
    """
    os.makedirs("cut", exist_ok=True)
    successful_cuts = 0

    for i in range(1, num_frames + 1):
        dump_file = f"{file_prefix}{i:04d}.dump"
        if preprocess_single_frame(dump_file, radius, cut_molecules):
            successful_cuts += 1

    print(f"Processing complete: {successful_cuts}/{num_frames} cuts performed.")


def train_potential(retrain=False):
    """Starts potential training. If retrain=True, training starts from scratch."""
    if retrain:
        os.system(
            f"sbatch -J train -o logs/train.out -e logs/train.err -N 1 --exclusive scripts/retrain.sh --time=3:00:00"
        )
    else:
        os.system(
            f"sbatch -J train -o logs/train.out -e logs/train.err -N 1 --exclusive scripts/train.sh --time=3:00:00"
        )
    while not Path("results/is_training_finished.txt").exists():
        time.sleep(5)


def process_frames(data_folder, processed_folder, num_frames, radius, cut):
    """Splits and processes frames from dump files."""
    dump_files = sorted(Path(data_folder).glob("*.dump"))
    first = True
    for dump_file in dump_files:
        print(f"Processing trajectory file: {dump_file}, splitting...")
        cut_dump(dump_file, f"{processed_folder}/frame_")
        if first:
            preprocess_multiple_frames(
                num_frames, f"{processed_folder}/frame_", radius, cut
            )
            first = False

    all_frames = sorted(Path(processed_folder).glob("frame_*.dump"))
    for frame in all_frames:
        convert("dump", "cfg", frame, frame.with_suffix(".cfg"), None)

    with open("results/candidates.cfg", "w") as candidates_file:
        for frame_cfg in sorted(Path(processed_folder).glob("frame_*.cfg")):
            with open(frame_cfg, "r") as f:
                candidates_file.write(f.read())


def run_calculations(frame_folder):
    """Generates input data and runs QM/MM calculations."""
    n_calcs = 0
    os.makedirs("final_cfgs", exist_ok=True)
    final_cfg_folder = Path("final_cfgs")
    for xyz_file in sorted(Path(frame_folder).glob("*cut*.xyz")):
        frame_id = xyz_file.stem.split("_")[-1]
        frame_qm_folder = Path(f"qm_calculations/frame_{frame_id}")
        frame_mm_folder = Path(f"mm_calculations/frame_{frame_id}")
        frame_qm_folder.mkdir(parents=True, exist_ok=True)
        frame_mm_folder.mkdir(parents=True, exist_ok=True)

        dump_file = xyz_file.with_suffix(".dump")

        generate_orca_input(
            "templates/orca_template.inp",
            xyz_file,
            frame_qm_folder / xyz_file.with_suffix(".inp"),
        )
        shutil.copy2(xyz_file, frame_qm_folder)

        cut_data(
            "templates/system.data",
            dump_file,
            frame_mm_folder / f"{xyz_file.name}.data",
        )

        os.system(
            f"cd {frame_qm_folder} && sbatch ../../scripts/run_qm.sh {xyz_file.name}.inp"
        )
        os.system(
            f"cd {frame_mm_folder} && sbatch ../../scripts/run_mm.sh {xyz_file.name}.data"
        )

    while True:
        n_qm_finished = len(Path("qm_calculations").glob("/*/is_qm_finished.txt"))
        n_mm_finished = len(Path("mm_calculations").glob("/*/is_mm_finished.txt"))

        if (n_qm_finished == n_calcs) and (n_mm_finished == n_calcs):
            break
        else:
            time.sleep(10)

    for dump_file in sorted(Path("mm_calculations").glob("/*/*.dump")):
        shutil.copy2(dump_file, final_cfg_folder / dump_file.name)

    for log_file in sorted(Path("qm_calculations").glob("/*/*.log")):
        shutil.copy2(log_file, final_cfg_folder / log_file.name)

    for log_file in sorted(final_cfg_folder.glob("*.log")):
        convert("orcadump", "cfg", log_file, log_file.with_suffix("cfg"), None)

    with open("results/train.cfg", "a") as train_file:
        for cfg_file in sorted(final_cfg_folder.glob("*.cfg")):
            with open(cfg_file, "r") as f:
                train_file.write(f.read())


def cut_and_calculate(frame_folder, radius, cutoff):
    """Cuts environment and performs additional QM/MM calculations."""
    for cfg_file in sorted(Path(frame_folder).glob("frame_*_cut.cfg")):
        config = Configuration()
        config.read_cfg(cfg_file, store_grade=True)
        max_atom = max(config.atom_data, key=lambda atom: atom["nbh_grades"])
        atom_coords = [max_atom["cartes_x"], max_atom["cartes_y"], max_atom["cartes_z"]]
        cut_nbh(
            cfg_file, cfg_file.with_suffix("_cut_max.cfg"), atom_coords, radius, cutoff
        )

    run_calculations(frame_folder, None)  # Run on all available frames


def process_and_train(
    data_folder,
    processed_folder,
    num_frames=100,
    radius=5.0,
    cut=True,
    grade_threshold=2.0,
):
    """Main process: frame processing, training, filtering, and iterations."""
    process_frames(data_folder, processed_folder, num_frames, radius, cut)
    run_calculations(processed_folder)
    return

    while True:
        os.system(
            "sbatch -J grade -o logs/grade.out -e logs/grade.err -N 1 -n 16 scripts/grade.sh"
        )
        while not Path("results/selected.cfg").exists():
            time.sleep(5)

        filter_by_grade("results/train.cfg", "results/filtered.cfg", grade_threshold)

        if (
            not Path("results/filtered.cfg").exists()
            or os.stat("results/filtered.cfg").st_size == 0
        ):
            break

        os.system(
            "sbatch -J selection -o logs/selection.out -e logs/selection.err -N 1 -n 8 scripts/select.sh"
        )
        while not Path("results/is_selection_finished.txt").exists():
            time.sleep(5)

        cut_and_calculate(processed_folder, radius, cut)

        retrain = check_cfg_ratio(
            "results/new_selected.cfg", "results/train.cfg", fraction=0.5
        )
        os.system("cat results/new_selected.cfg >> results/train.cfg")
        print(
            f"Updating train.cfg: {'Retraining' if retrain else 'Fine-tuning'} the potential"
        )
        train_potential(retrain=retrain)

    print("Process completed")


if __name__ == "__main__":
    process_and_train(
        "data",
        "processed_frames",
        num_frames=100,
        radius=5.0,
        cut=True,
        grade_threshold=2.0,
    )


"""
project_root/
│── data/                # Raw input dump files
│   ├── 1.dump
│   ├── 2.dump
│   └── ...
│── processed_frames/    # Stores processed frames from dumps
│   ├── frame_0001.dump
│   ├── frame_0002.dump
│   ├── ...
|   └── cut/
│       ├── frame_0001_cut.dump
│       ├── frame_0001_cut.xyz
│       └── ...
│── qm_calculations/     # QM calculation folders
│   ├── frame_0001/
│   ├── frame_0002/
│   └── ...
│── mm_calculations/     # MM calculation folders
│   ├── frame_0001/
│   ├── frame_0002/
│   └── ...
│── final_cfgs/          # Stores final configuration files
│── scripts/             # Stores executable scripts (e.g., ORCA, MM, grading, training)
│   ├── run_orca.sh
│   ├── run_mm.sh
│   ├── grade.sh
│   ├── train.sh
│   ├── select.sh
│── templates/           # Stores reusable templates (e.g., ORCA input files)
│   ├── orca_template.inp
│── results/             # Stores training results
│   ├── train.cfg
│   ├── filtered.cfg
│   ├── selected.cfg
│── logs/                # Stores output logs
│   ├── train.out
│   ├── selection.out
│   └── ...
└── main.py              # Main script (this script)

"""
