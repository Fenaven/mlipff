import os
import time
from pathlib import Path
from utils import filter_by_grade
from frames_processing import process_multiple_frames


def train_potential():
    """Запускает обучение потенциала."""
    os.system(
        "sbatch -J train -o train.out -e train.err -N 1 --exclusive ./train.sh --time=3:00:00"
    )
    while not Path("is_training_finished.txt").exists():
        time.sleep(5)


def process_and_train(
    md_folder, num_frames=100, radius=5.0, cutoff=True, grade_threshold=2.0
):
    """
    Основной процесс: обработка кадров, обучение, фильтрация, итерации.
    """
    process_multiple_frames(num_frames, md_folder + "/frame_", radius, cutoff)

    with open("train.cfg", "w") as train_file:
        for i in range(1, num_frames + 1):
            cfg_file = f"{md_folder}/frame_{i:04d}_cut.cfg"
            if Path(cfg_file).exists():
                with open(cfg_file, "r") as f:
                    train_file.write(f.read())

    train_potential()

    while True:
        filter_by_grade("train.cfg", "filtered.cfg", grade_threshold)

        if not Path("filtered.cfg").exists() or os.stat("filtered.cfg").st_size == 0:
            break

        os.system(
            "sbatch -J selection -o selection.out -e selection.err -N 1 -n 8 ./select.sh"
        )
        while not Path("is_selection_finished.txt").exists():
            time.sleep(5)

        os.rename("filtered.cfg", "train.cfg")
        train_potential()

    print("Процесс завершён")


if __name__ == "__main__":
    process_and_train(
        "MD_dumps", num_frames=100, radius=5.0, cutoff=True, grade_threshold=2.0
    )
