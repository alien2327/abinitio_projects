"""
Reference: https://legacy.materialsproject.org/materials/mp-18725/

Final Magnetic Moment: 0.000 μB
Magnetic Ordering: AFM
Formation Energy / Atom: -3.212 eV
Energy Above Hull / Atom: 0.000 eV
Density: 5.48 g/cm3
Decomposes To: Stable
Band Gap: 2.338 eV

Run type: GGA+U
Cutoff E: 700 eV
k-points: None
U values: Cr: 3.7 eV
Final E/Atom: -8.5802 eV
Corrected E: -187.8443 eV
    - Uncorrected energy = -171.6043 eV
    - Composition-based energy adjustment (-0.687 eV/atom x 12.0 atoms) = -8.2440 eV
        - MP2020 anion correction (oxide)
    - Composition-based energy adjustment (-1.999 eV/atom x 4.0 atoms) = -7.9960 eV
        - MP2020 GGA/GGA+U mixing correction (Cr)
"""
import warnings
warnings.filterwarnings('ignore')

import sys
import numpy as np

from ase import units
from ase.io import write
from ase import Atom, Atoms

from gpaw import GPAW, PW, FermiDirac

temp = int(sys.argv[1])

a, b, c = 5.317, 5.615, 7.667
a_t, b_t, c_t = 90, 90, 90

Y0 = Atom('Y', (a*0.0189, b*0.9326, c*0.7500))
Y1 = Atom('Y', (a*0.4811, b*0.4326, c*0.7500))
Y2 = Atom('Y', (a*0.5189, b*0.5674, c*0.2500))
Y3 = Atom('Y', (a*0.9811, b*0.0674, c*0.2500))

Cr0 = Atom('Cr', (a*0.0000, b*0.5000, c*0.0000))
Cr1 = Atom('Cr', (a*0.0000, b*0.5000, c*0.5000))
Cr2 = Atom('Cr', (a*0.5000, b*0.0000, c*0.5000))
Cr3 = Atom('Cr', (a*0.5000, b*1.0000, c*0.0000))

O0 =  Atom('O', (a*0.1119, b*0.4591, c*0.2500))
O1 =  Atom('O', (a*0.1930, b*0.1976, c*0.5574))
O2 =  Atom('O', (a*0.1930, b*0.1976, c*0.9426))
O3 =  Atom('O', (a*0.3070, b*0.6976, c*0.9426))
O4 =  Atom('O', (a*0.3070, b*0.6976, c*0.5574))
O5 =  Atom('O', (a*0.3881, b*0.9591, c*0.2500))
O6 =  Atom('O', (a*0.6119, b*0.0409, c*0.7500))
O7 =  Atom('O', (a*0.6930, b*0.3024, c*0.4426))
O8 =  Atom('O', (a*0.6930, b*0.3024, c*0.0574))
O9 =  Atom('O', (a*0.8070, b*0.8024, c*0.0574))
O10 = Atom('O', (a*0.8070, b*0.8024, c*0.4426))
O11 = Atom('O', (a*0.8881, b*0.5409, c*0.7500))

YCrO3 = Atoms(
    symbols=[
        Y0, Y1, Y2, Y3,
        Cr0, Cr1, Cr2, Cr3,
        O0, O1, O2, O3, O4, O5, O6, O7, O8, O9, O10, O11
    ],
    cell=[
        a, b, c,
        a_t, b_t, c_t
    ],
    pbc=True
)

YCrO3.calc = GPAW(
    mode='lcao',
    basis='dzp',
    xc='LDA',
    hund=True,
    spinpol=True,
    occupations={
        'name': 'fermi-dirac',
        'width': temp * units.kB
    },
    setups={
        'Cr': ':d,3.7'
    },
    parallel = {
        'augment_grids': True,
        'sl_auto': True,
        'use_elpa': True
    },
    maxiter=1000,
    kpts=(3, 3, 3),
    convergence={
        # 'density': 2.0e-4,
        # 'eigenstates': 2.0e-7,
        'bands': 'occupied'
    },
    txt=f'log/YCrO3_{temp}K_scf.log'
)

YCrO3.get_potential_energy()
YCrO3.calc.write(f'log/YCrO3_{temp}K_scf.gpw', mode='all')
