import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.io import read
from ase.units import GPa, kg
from calorine.calculators import CPUNEP
from calorine.tools import get_elastic_stiffness_tensor, get_force_constants, relax_structure
from matplotlib import pyplot as plt
from pandas import DataFrame


structure = read('model.xyz')
calculator = CPUNEP('nepfile')
structure.calc = calculator
relax_structure(structure, fmax=0.0001)

cij_rlx = get_elastic_stiffness_tensor(structure, fmax=1e-5)
cij_clamped = get_elastic_stiffness_tensor(structure, clamped=True)

with np.printoptions(precision=1, suppress=True):
    print('Relaxed elastic constants')
    print(cij_rlx)
    print('')
    print('Clamped elastic constants')
    print(cij_clamped)
