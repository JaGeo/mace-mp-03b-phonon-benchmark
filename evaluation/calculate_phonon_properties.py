import os
import shutil
import sys

sys.path.append("functions")
from functions import fine_tuning
from functions.phonon import Phonon_Properties
from functions.plots import PhononPlotter
import numpy as np
import matplotlib.pyplot as splt
import pandas as pd

for mp_id in os.listdir("./dft_pbe_phonons"):
    try:
        df = pd.read_csv("./phonon_properties.csv")
        if mp_id in df["material"].values:
            continue
    except:
        pass
    print(mp_id)
    phonopy_yaml = (
        f"./mp_0b3_medium_phonons/{mp_id}/phonopy.yaml"
    )
    phonopy_yaml_dft = (
        f"./dft_pbe_phonons/{mp_id}/phonopy.yaml"
    )

    if not os.path.exists(phonopy_yaml) or not os.path.exists(phonopy_yaml_dft):
        print(f"phonopy.yaml not found for {mp_id}")
        continue

    phonon_prop = Phonon_Properties(phonopy_yaml=phonopy_yaml)
    if os.path.exists(
        os.path.join("./dft_pbe_phonons", mp_id, "FORCE_SETS")
    ):
        phonon_prop_dft = Phonon_Properties(
            phonopy_yaml=phonopy_yaml_dft,
            force_sets_filename=os.path.join(
                "./dft_pbe_phonons", mp_id, "FORCE_SETS"
            ),
        )
    else:
        phonon_prop_dft = Phonon_Properties(phonopy_yaml=phonopy_yaml_dft)

    try:
        phonon_distance = phonon_prop._run_bands_structure_dict()["distances"]
        phonon_distance_dft = phonon_prop_dft._run_bands_structure_dict()["distances"]
        phonon_frequencies = phonon_prop._run_bands_structure_dict()["frequencies"]
        phonon_frequencies_dft = phonon_prop_dft._run_bands_structure_dict()[
            "frequencies"
        ]
        labels, connections = phonon_prop_dft._get_label_and_connection()
        plot = PhononPlotter(
            distances_set=([phonon_distance_dft, phonon_distance_dft]),
            frequencies_set=([phonon_frequencies_dft, phonon_frequencies]),
            x_labels=labels,
            connections=connections,
            colors=["black", "red"],
            legend_labels=["DFT PBE", "MP-0b3"],
        )
        fig, ax = plot.beautiful_phonon_plotter()
        plt.savefig(
            f"./mace_phonons/plots/{mp_id}.png",
            dpi=300,
            bbox_inches="tight",
        )
    except:
        print(f"PLOT NOT POSSIBLE FOR {mp_id}")

    rmse, mae, r2, phon_rmse, phon_mae, phon_scaled_rmse, phon_rrmse = (
        phonon_prop.get_metrics(other=phonon_prop_dft)
    )
    min_f = phonon_prop.get_min_frequency()
    max_f = phonon_prop.get_max_frequency()
    mean_f = phonon_prop.get_mean_frequency()
    volume = phonon_prop.get_volume()
    rmsd, maxd = phonon_prop.get_displacements(other=phonon_prop_dft)
    spacegroup = phonon_prop.get_space_group()

    data = {
        "material": mp_id,
        "rmse": rmse,
        "mae": mae,
        "r2": r2,
        "phonon_rmse": phon_rmse,
        "phonon_mae": phon_mae,
        "phonon_scaled_rmse": phon_scaled_rmse,
        "phonon_rrmse": phon_rrmse,
        "rmsd": rmsd,
        "maxd": maxd,
        "spacegroup": spacegroup,
        "min_frequency": min_f,
        "max_frequency": max_f,
        "mean_frequency": mean_f,
        "volume": volume,
        "benchmark_spacegroup": phonon_prop_dft.get_space_group(),
        "benchmark_min_frequency": phonon_prop_dft.get_min_frequency(),
        "benchmark_max_frequency": phonon_prop_dft.get_max_frequency(),
        "benchmark_mean_frequency": phonon_prop_dft.get_mean_frequency(),
        "benchmark_volume": phonon_prop_dft.get_volume(),
    }

    if os.path.exists("./phonon_properties.csv"):
        df = pd.concat([df, pd.DataFrame(data, index=[0])], ignore_index=True)
    else:
        df = pd.DataFrame(data, index=[0])

    print(df)
    df.to_csv("./phonon_properties.csv", index=False)
