import warnings
warnings.filterwarnings('ignore')

import sys
import numpy as np

from gpaw import restart

temp = int(sys.argv[1])
index = int(sys.argv[2])

a, b, c = 5.317, 5.615, 7.667
a_t, b_t, c_t = 90, 90, 90

YCaCr2O3, calc = restart(f'log/YCaCr2O3_{temp}K_{index}_scf.gpw')

bandpath = YCaCr2O3.cell.bandpath(npoints=512)
YCaCr2O3.calc = YCaCr2O3.calc.fixed_density(
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
    txt=f'log/YCaCr2O3_{temp}K_{index}_nscf.log'
)

YCaCr2O3.calc.band_structure()
YCaCr2O3.calc.get_dos(spin=0, npts=2001)
YCaCr2O3.calc.get_dos(spin=1, npts=2001)
YCaCr2O3.calc.write(f'log/YCaCr2O3_{temp}K_{index}_nscf.gpw', mode='all')
