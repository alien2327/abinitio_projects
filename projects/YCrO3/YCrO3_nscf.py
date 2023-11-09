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

from gpaw import restart

temp = int(sys.argv[1])

YCrO3, calc = restart(f'log/YCrO3_{temp}K_scf.gpw')

bandpath = YCrO3.cell.bandpath(npoints=512)
YCrO3.calc = YCrO3.calc.fixed_density(
    nbands=120,
    symmetry='off',
    kpts={
        'kpts': bandpath.kpts,
        'path': bandpath.path,
        'npoints': 60
    },
    convergence={
        # 'density': 2.0e-4,
        # 'eigenstates': 2.0e-7,
        'bands': 'all'
    },
    txt=f'log/YCrO3_{temp}K_nscf.log'
)

YCrO3.calc.band_structure()
YCrO3.calc.get_dos(spin=0, npts=2001)
YCrO3.calc.get_dos(spin=1, npts=2001)
YCrO3.calc.write(f'log/YCrO3_{temp}K_nscf.gpw', mode='all')
