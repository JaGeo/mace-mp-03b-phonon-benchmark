import numpy as np
import os
import matplotlib.pyplot as plt
from ase.io import read, write
from sklearn.metrics import r2_score, root_mean_squared_error, mean_absolute_error


def calculate_rmse(true_value, evaluated_value):
    """calculate root mean square error
    true_value: np.ndarray - true value
    evalutated_value: np.ndarray - evaluated value
    """
    if len(true_value) != len(evaluated_value):
        return None
    elif len(true_value) == 0:
        raise ValueError(f"Length of true value and evaluated value is 0")
    else:
        diff = np.array(true_value) - np.array(evaluated_value)
        rmse = np.sqrt(np.mean(diff**2))
        return rmse


def calculate_rrmse(true_value, evaluated_value):
    """calculate relative root mean square error
    true_value: np.ndarray - true value
    evaluated_value: np.ndarray - evaluated value
    """
    if len(true_value) != len(evaluated_value):
        return None
    elif len(true_value) == 0:
        raise ValueError(f"Length of true value and evaluated value is 0")
    else:
        rms = np.sqrt(np.mean(true_value**2))
        rrmse = calculate_rmse(true_value, evaluated_value) / rms
        return rrmse


def calculate_mae(true_value, evaluated_value):
    """calculate mean absolute error
    true_value: np.ndarray - true value
    evaluated_value: np.ndarray - evaluated value
    """
    if len(true_value) != len(evaluated_value):
        return None
    elif len(true_value) == 0:
        raise ValueError(f"Length of true value and evaluated value is 0")
    else:
        diff = np.array(true_value) - np.array(evaluated_value)
        mae = np.mean(np.abs(diff))
        return mae


def calculate_r2(true_value, evaluated_value):
    """calculate r2 score
    true_value: np.ndarray - true value
    evaluated_value: np.ndarray - evaluated value
    """
    if len(true_value) != len(evaluated_value):
        return None
    elif len(true_value) == 0:
        raise ValueError(f"Length of true value and evaluated value is 0")
    else:
        r2 = r2_score(
            np.array(true_value).flatten(), np.array(evaluated_value).flatten()
        )
        return r2


def phonon_rmse(true_mesh_dict, evaluated_mesh_dict):
    """calculate root mean square error especially for phonons
    over the whole mesh in the Brillouin zone
    true_mesh_dict: dict - true mesh dict  phonon.get_mesh_dict()
    evaluated_mesh_dict: dict - evaluated mesh dict phonon.get_mesh_dict()
    """
    f_dft = true_mesh_dict["frequencies"]
    f_ml = evaluated_mesh_dict["frequencies"]
    w = true_mesh_dict["weights"]
    f_dft_sorted = np.sort(f_dft, axis=1)
    f_ml_sorted = np.sort(f_ml, axis=1)
    w_reshape = np.repeat(w[:, np.newaxis], len(f_dft_sorted[0, :]), axis=1)
    diff = f_dft_sorted - f_ml_sorted
    # rmse = np.sqrt(np.mean(w_reshape*diff**2)/np.mean(w_reshape))
    # same expression but more explicit
    rmse = np.sqrt(1 / len(f_dft_sorted[0]) * (np.sum(w_reshape * diff**2) / np.sum(w)))
    return rmse


def phonon_mae(true_mesh_dict, evaluated_mesh_dict):
    """calculate mean absolute error especially for phonons
    over the whole mesh in the Brillouin zone
    true_mesh_dict: dict - true mesh dict  phonon.get_mesh_dict()
    evaluated_mesh_dict: dict - evaluated mesh dict phonon.get_mesh_dict()
    """
    f_dft = true_mesh_dict["frequencies"]
    f_ml = evaluated_mesh_dict["frequencies"]
    w = true_mesh_dict["weights"]
    f_dft_sorted = np.sort(f_dft, axis=1)
    f_ml_sorted = np.sort(f_ml, axis=1)
    w_reshape = np.repeat(w[:, np.newaxis], len(f_dft_sorted[0, :]), axis=1)
    diff = np.abs(f_dft_sorted - f_ml_sorted)
    mae = 1 / len(f_dft_sorted[0]) * np.sum(w_reshape * diff) / np.sum(w)
    return mae


def phonon_rrmse(true_mesh_dict, evaluated_mesh_dict):
    """calculate relative root mean square error especially for phonons
    over the whole mesh in the Brillouin zone
    true_mesh_dict: dict - true mesh dict  phonon.get_mesh_dict()
    evaluated_mesh_dict: dict - evaluated mesh dict phonon.get_mesh_dict()
    """
    f_dft = true_mesh_dict["frequencies"]
    w = true_mesh_dict["weights"]
    f_dft_sorted = np.sort(f_dft, axis=1)
    w_reshape = np.repeat(w[:, np.newaxis], len(f_dft_sorted[0, :]), axis=1)
    f_dft_sorted = np.sort(f_dft, axis=1)
    rmse = phonon_rmse(true_mesh_dict, evaluated_mesh_dict)
    rms = np.sqrt(
        1 / len(f_dft_sorted[0]) * (np.sum(w_reshape * f_dft_sorted**2) / np.sum(w))
    )
    phonon_rrmse = rmse / rms
    return phonon_rrmse


def phonon_shifted_rmse(true_mesh_dict, evaluated_mesh_dict):
    """calculate root mean square error especially for phonons
    over the whole mesh in the Brillouin zone
    true_mesh_dict: dict - true mesh dict  phonon.get_mesh_dict()
    evaluated_mesh_dict: dict - evaluated mesh dict phonon.get_mesh_dict()
    """
    f_dft = true_mesh_dict["frequencies"]
    f_ml = evaluated_mesh_dict["frequencies"]
    scaling_factor = np.max(f_dft) / np.max(f_ml)
    f_ml = f_ml * scaling_factor
    evaluated_mesh_dict["frequencies"] = f_ml
    rmse = phonon_rmse(true_mesh_dict, evaluated_mesh_dict)
    return rmse


def mean_phonon_frequency(mesh_dict):
    """calculate mean phonon frequency
    mesh_dict: dict - mesh dict  phonon.get_mesh_dict()
    """
    f = mesh_dict["frequencies"]  # 2D-Array (N_q x N_modes)
    w = mesh_dict["weights"]  # 1D-Array (N_q,)
    w_expanded = np.repeat(w[:, np.newaxis], f.shape[1], axis=1)
    mask = f > 0
    f_filtered = f[mask]
    w_filtered = w_expanded[mask]
    f_mean = np.sum(f_filtered * w_filtered) / np.sum(w_filtered)
    return f_mean
