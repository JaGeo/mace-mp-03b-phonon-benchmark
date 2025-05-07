from pymatgen.io.vasp import inputs
from mp_api.client import MPRester
from atomate2.vasp.flows.phonons import PhononMaker
from atomate2.vasp.jobs.mp import MPGGAStaticMaker
from atomate2.vasp.powerups import (
    update_user_incar_settings,
    update_user_kpoints_settings,
)
from atomate2.vasp.flows.mp import MPGGADoubleRelaxMaker
from pymatgen.io.vasp import Incar, Kpoints
from atomate2.common.jobs.phonons import get_supercell_size
from atomate2.common.flows.phonons import BasePhononMaker
from atomate2.forcefields.jobs import (
    MACERelaxMaker,
    MACEStaticMaker,
    ForceFieldRelaxMaker,
    ForceFieldStaticMaker,
)
from jobflow import run_locally
import os
from pymatgen.core import Structure
import json
from phonopy import load
from phonopy.units import VaspToTHz
from pymatgen.io.phonopy import get_pmg_structure, get_ph_bs_symm_line
import numpy as np
from pymatgen.symmetry.kpath import KPathSeek
import copy
from phonopy.phonon.band_structure import get_band_qpoints
from phonopy.phonon.band_structure import BandStructure, get_band_qpoints_by_seekpath
import sys
from functions.metric import (
    calculate_rmse,
    calculate_mae,
    calculate_r2,
    phonon_rmse,
    phonon_rrmse,
    phonon_mae,
    phonon_shifted_rmse,
    mean_phonon_frequency,
)
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
import json
import os
from pymatgen.analysis.structure_matcher import StructureMatcher


def calculate_phonon(
    model_path: str,
    structure: Structure,
    run_dir: str,
    name: str,
    supercell_matrix: None,
    supercell_size: int = 20,
    symprec: float = 1e-4,
    displacement: float = 0.01,
    device: str = "cuda",
):
    """calculate phonon dispersion relation with phonopy and mace model"""

    calculator_kwargs = {
        "device": "cpu",
        "default_dtype": "float64",
        "model": model_path,
    }

    relax_kwargs = {
        "fmax": 0.00001,
        "steps": 2000,
    }

    phonon = PhononMaker(
        name=name,
        bulk_relax_maker=ForceFieldRelaxMaker(
            calculator_kwargs=calculator_kwargs, force_field_name="MACE"
        ),
        born_maker=None,
        phonon_displacement_maker=ForceFieldStaticMaker(
            calculator_kwargs=calculator_kwargs, force_field_name="MACE"
        ),
        static_energy_maker=ForceFieldStaticMaker(
            calculator_kwargs=calculator_kwargs, force_field_name="MACE"
        ),
        min_length=supercell_size,
        displacement=displacement,
        symprec=symprec,
    ).make(structure, supercell_matrix=supercell_matrix)
    run_locally(phonon, root_dir=run_dir, create_folders=False)


class Phonon_Properties:
    def __init__(
        self, phonopy_yaml: str, force_sets: str = None, symprec: float = 1e-4
    ):

        self.phonopy_yaml = phonopy_yaml
        self.force_sets = force_sets
        self.symprec = symprec

        # Initialize missing attributes
        self.bandstructure_dict = None
        self.dos = None
        self.phonon_object = None
        self.structure = None
        self.kpath = None
        self.path = None
        self.band_qpoints = None
        self.labels_for_plot = None
        self.connections = None
        self.space_group = None
        self.mesh_dict = None
        self.dos_dict = None

    def _load_phonopy_object(self):
        if self.phonon_object is None:
            if self.force_sets is None:
                self.phonon_object = load(
                    phonopy_yaml=self.phonopy_yaml,
                    factor=VaspToTHz,
                    is_nac=False,
                    symprec=self.symprec,
                )
            else:
                self.phonon_object = load(
                    phonopy_yaml=self.phonopy_yaml,
                    factor=VaspToTHz,
                    is_nac=False,
                    symprec=self.symprec,
                    force_sets_filename=self.force_sets,
                )
            self.phonon_object._run_force_constants_from_forces()
            fc = self.phonon_object.force_constants
            self.phonon_object.set_force_constants(force_constants=fc)
            self.phonon_object._set_dynamical_matrix()
        return self.phonon_object

    def _get_structure(self):
        if self.structure is None:
            phonon = self._load_phonopy_object()
            self.structure = get_pmg_structure(phonon.primitive)
        return self.structure

    def _compute_kpath(self):
        """kpath calculated with seekpath"""
        if self.path or self.kpath is None:
            structure = self._get_structure()
            highsymmkpath = KPathSeek(structure=structure, symprec=self.symprec)
            self.kpath = highsymmkpath._kpath
            self.path = copy.deepcopy(self.kpath["path"])
            for idx, labelset in enumerate(self.kpath["path"]):
                for i, label in enumerate(labelset):
                    self.path[idx][i] = self.kpath["kpoints"][label]
        return self.kpath["kpoints"], self.path

    def _get_band_qpoints(self):
        # Berechne die Bandpunkte
        if self.band_qpoints is None:
            band_paths = self._compute_kpath()[1]
            self.band_qpoints = get_band_qpoints(band_paths=band_paths, npoints=101)
        return self.band_qpoints

    def _get_label_and_connection(self):
        if self.labels_for_plot or self.connections is None:
            structure = self._load_phonopy_object().primitive
            npoints = 101
            self.bands, self.labels_for_plot, self.connections = (
                get_band_qpoints_by_seekpath(
                    primitive=structure, npoints=npoints, is_const_interval=True
                )
            )
        return self.labels_for_plot, self.connections

    def _compute_properties(self):
        # Errechne den k-Pfad und erhalte ein Dictionary und den konkreten Pfad
        self.kpath_dict, self.kpath_concrete = self._compute_kpath(
            structure=self._get_structure(),
            kpath_scheme=self.kpath_scheme,
            symprec=self.symprec,
        )

    def _run_bands_structure_dict(self, paths=None, labels=None):
        if self.bandstructure_dict is None:
            phonon = self._load_phonopy_object()
            if paths is None:
                paths = self._get_band_qpoints()
            if labels is None:
                labels = self._get_label_and_connection()[0]
            phonon.run_band_structure(
                paths=paths, labels=labels, is_band_connection=True
            )
            self.labels_for_plot = labels
            self.paths = paths
            self.bandstructure_dict = phonon.get_band_structure_dict()
        return self.bandstructure_dict

    def _run_dos(self, mesh=None, symmetry=True):
        phonon = self._load_phonopy_object()
        if mesh is None:
            mesh = [20, 20, 20]
        phonon.run_mesh(mesh=mesh, is_mesh_symmetry=symmetry, is_time_reversal=symmetry)
        phonon.run_total_dos(sigma=0.05)
        self.dos_dict = phonon.get_total_dos_dict()
        self.mesh_dict = phonon.get_mesh_dict()
        return self.dos_dict

    def get_space_group(self):
        if self.space_group is None:
            phonon = self._load_phonopy_object()
            self.space_group = phonon.symmetry.get_international_table()
        return self.space_group

    def _calculate_metrics(self, other):
        if self.bandstructure_dict is None:
            self.bandstructure_dict = self._run_bands_structure_dict()
        if other.bandstructure_dict is None:
            other.bandstructure_dict = other._run_bands_structure_dict()
        freq1 = np.array(self.bandstructure_dict["frequencies"])
        freq2 = np.array(other.bandstructure_dict["frequencies"])
        self.rmse = calculate_rmse(true_value=freq1, evaluated_value=freq2)
        self.mae = calculate_mae(true_value=freq1, evaluated_value=freq2)
        self.r2 = calculate_r2(true_value=freq1, evaluated_value=freq2)

        if self.dos_dict is None:
            self.dos_dict = self._run_dos()
        if other.dos_dict is None:
            other.dos_dict = other._run_dos()
        mesh1 = self.mesh_dict
        mesh2 = other.mesh_dict
        if mesh1["frequencies"].shape != mesh2["frequencies"].shape:
            self._run_dos(symmetry=False)
            other._run_dos(symmetry=False)
        mesh1 = self.mesh_dict
        mesh2 = other.mesh_dict
        self.phonon_rmse = phonon_rmse(true_mesh_dict=mesh1, evaluated_mesh_dict=mesh2)
        self.phonon_mae = phonon_mae(true_mesh_dict=mesh1, evaluated_mesh_dict=mesh2)
        self.phonon_scaled_rmse = phonon_shifted_rmse(
            true_mesh_dict=mesh1, evaluated_mesh_dict=mesh2
        )
        self.phonon_rrmse = phonon_rrmse(
            true_mesh_dict=mesh1, evaluated_mesh_dict=mesh2
        )

    def get_metrics(self, other):
        """calculate metrics for phonon dispersion relation and dos
        other: Phonon_Properties - other phonon properties object
        """
        self._calculate_metrics(other)
        return (
            self.rmse,
            self.mae,
            self.r2,
            self.phonon_rmse,
            self.phonon_mae,
            self.phonon_scaled_rmse,
            self.phonon_rrmse,
        )

    def get_max_frequency(self):
        """get max frequency from phonon dispersion relation"""
        if self.bandstructure_dict is None:
            self.bandstructure_dict = self._run_bands_structure_dict()
        self.max = np.max(np.array(self.bandstructure_dict["frequencies"]))
        return self.max

    def get_min_frequency(self):
        """get min frequency from phonon dispersion relation"""
        if self.bandstructure_dict is None:
            self.bandstructure_dict = self._run_bands_structure_dict()
        self.min = np.min(np.array(self.bandstructure_dict["frequencies"]))
        return self.min

    def get_mean_frequency(self):
        """get mean frequency from phonon dispersion relation"""
        if self.dos_dict is None:
            self.dos_dict = self._run_dos()
        mesh = self.mesh_dict
        self.mean = mean_phonon_frequency(mesh_dict=mesh)
        return self.mean

    def get_volume(self):
        """get volume of the structure"""
        if self.structure is None:
            self.structure = self._get_structure()
        self.volume = self.structure._lattice.volume
        return self.volume

    def get_displacements(self, other):
        """get structure properties from phonon dispersion relation"""
        if self.structure is None:
            self.structure = self._get_structure()
        if other.structure is None:
            other.structure = other._get_structure()
        try:
            self.rms_displacement, self.max_displacement = (
                StructureMatcher().get_rms_dist(self.structure, other.structure)
            )
        except:
            self.rms_displacement = None
            self.max_displacement = None
        return self.rms_displacement, self.max_displacement
