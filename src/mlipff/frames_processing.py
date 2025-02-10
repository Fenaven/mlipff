import random
from utils import cut_nbh
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


def process_single_frame(dump_file, radius, cut):
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

    output_name = dump_file.replace(".dump", "_cut")
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


def process_multiple_frames(num_frames, file_prefix, radius, cut_molecules):
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
    successful_cuts = 0

    for i in range(1, num_frames + 1):
        dump_file = f"{file_prefix}{i:04d}.dump"
        if process_single_frame(dump_file, radius, cut_molecules):
            successful_cuts += 1

    print(f"Processing complete: {successful_cuts}/{num_frames} cuts performed.")
