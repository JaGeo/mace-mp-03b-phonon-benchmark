from jobflow import run_locally
from pymatgen.core import Structure
from pymatgen.transformations.advanced_transformations import CubicSupercellTransformation
from atomate2.common.jobs.phonons import generate_phonon_displacements, get_supercell_size
from mp_api.client import MPRester


#input POSCAR Structure
structure = MPRester("fkFg0Vn7Al6jNXPvNdcWdBaiQTeqnL2u").get_structure_by_material_id('mp-2809').to_conventional()
structure.to(filename='POSCAR',fmt='poscar')
#generate Cubic Super Cell and get the transformation matrix
transformation = CubicSupercellTransformation(min_length=20)
supercell = transformation.apply_transformation(structure)
transformation_matrix = transformation.transformation_matrix.tolist()
supercell.to("SPOSCAR", fmt='poscar')
print(transformation_matrix)

#test mit phonon workflow code => Displacements mit diesem code nicht weniger, dafür ist die superzelle nicht quadratisch
#supercell_matrix = get_supercell_size(structure,min_length=20,prefer_90_degrees=True)
#transformation_mat = run_locally(supercell_matrix)

#displacements
displacements = generate_phonon_displacements(structure=structure, supercell_matrix=transformation_matrix, displacement=0.01, sym_reduce=True, symprec=1e-4, kpath_scheme="seekpath", code="vasp", use_symmetrized_structure=None)
responses = run_locally(displacements)

for inx, st in enumerate(responses[displacements.uuid][1].output):
    st.to(filename=f'POSCAR_{inx}',fmt='poscar')
