import sys

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import tqdm

from ase import units
from gpaw import restart

temp = int(sys.argv[1])

def get_data():
    YCrO3, calc = restart(f'log/YCrO3_{temp}K_nscf.gpw')

    atoms_indices = {}
    atoms_set = set(YCrO3.get_chemical_symbols())
    total_atoms = np.array(YCrO3.get_chemical_symbols())
    for atom in atoms_set:
        indices = np.where(total_atoms==atom)[0]
        atoms_indices[atom] = indices.tolist()
        
    pdos_0, pdos_1 = [], []

    pbar = tqdm.tqdm(total=len(YCrO3), desc=f'Total process')
    e_0, e_1 = None, None
    for atom, indicies in atoms_indices.items():
        pdos_a0_0, pdos_a0_1 = None, None
        
        for i, idx in enumerate(indicies):
            if i == 0:
                e_0, pdos_a0_0 = calc.get_orbital_ldos(spin=0, a=idx, npts=2001)
                e_1, pdos_a0_1 = calc.get_orbital_ldos(spin=1, a=idx, npts=2001)
            else:
                _, pdos_tmp_0 = calc.get_orbital_ldos(spin=0, a=idx, npts=2001)
                _, pdos_tmp_1 = calc.get_orbital_ldos(spin=1, a=idx, npts=2001)
                pdos_a0_0 += pdos_tmp_0
                pdos_a0_1 += pdos_tmp_1
            pbar.update(1)

        pdos_0.append(pdos_a0_0)
        pdos_1.append(pdos_a0_1)

    bs = calc.band_structure()
    dos_0 = np.sum(pdos_0, axis=0)
    dos_1 = np.sum(pdos_1, axis=0)
    e_f = calc.get_fermi_level()
    bs._energies -= e_f
    bs._reference = 0.0
    print("Fermi level: ", e_f)

    return atoms_indices, bs, e_0, e_1, e_f, dos_0, dos_1, pdos_0, pdos_1

def plot():
    atoms_indices, bs, e_0, e_1, e_f, dos_0, dos_1, pdos_0, pdos_1 = get_data()
    emin, emax = -5, 7

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.01
    rect_bs = [left, bottom, width, height]
    rect_dos = [left + width + spacing, bottom, 0.2, height]

    plt.figure(figsize=(12, 10))
    plt.rcParams.update({'font.size': 10})
    ax_bs = plt.axes(rect_bs)
    ax_bs.tick_params(direction='in', top=True)

    ax_dos = plt.axes(rect_dos)
    ax_dos.tick_params(direction='in', labelleft=False)

    ax_bs = bs.plot(
        ax=ax_bs,
        show=False,
        emin=emin,
        emax=emax
    )
    ax_bs.set_ylabel(r'$E - E_{\mathrm{Fermi}}$ (eV)')  
    ax_bs.set_xlabel(r'Wave vector')

    ax_dos.set_ylim(emin, emax)
    ax_dos.set_xlabel("DOS")

    ax_dos.plot(
        dos_0, e_0 - e_f,
        color='k',
        label="Total DOS (spin ↑)",
        zorder=1
    )
    ax_dos.plot(
        -dos_1, e_1 - e_f,
        color='k',
        linestyle='--',
        label="Total DOS (spin ↓)",
        zorder=2
    )
    ax_dos.fill_between(
        dos_0, e_0 - e_f,
        alpha=0.3,
        interpolate=True,
        color='grey',
        zorder=1,
    )
    ax_dos.fill_between(
        -dos_1, e_1 - e_f,
        alpha=0.3,
        interpolate=True,
        color='grey',
        zorder=2,
    )

    for idx, (atom, indicies) in enumerate(atoms_indices.items()):
        ax_dos.plot(pdos_0[idx], e_0 - e_f, label=f"{atom} PDOS spin ↑")
        ax_dos.plot(-pdos_1[idx], e_1 - e_f, label=f"{atom} PDOS spin ↓", linestyle='--')

    ax_dos.legend(
        loc='upper left', 
        bbox_to_anchor=(1.0, 0.5)
    )

    ax_dos.set_xlim(-25, 25)
    ax_dos.grid(
        color='silver',
        linestyle='--',
        linewidth=0.5,
        alpha=0.75,
        zorder=-1
    )

    plt.savefig(f'results/YCrO3_{temp}K_dos.png', bbox_inches='tight', dpi=300)

plot()