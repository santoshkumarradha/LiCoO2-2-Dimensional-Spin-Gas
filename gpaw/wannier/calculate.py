from ase import Atoms
from ase.build import bulk
from gpaw import GPAW, FermiDirac, PW
from ase.io import read
import os
import gpaw.wannier90 as w90
from gpaw import GPAW
from ase.parallel import parprint
import pickle
from datetime import datetime

parprint("Calculation starting {}".format(
    datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
if True:

    parprint('{} Reading cif \n'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    a = read("small.cif")

    parprint('{} Calculation starting \n'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    calc = GPAW(mode=PW(600),
                xc='LDA',
                occupations=FermiDirac(width=0.01),
                convergence={'density': 1.e-6},
                symmetry='off',
                nbands=-20,
                kpts={
                    'size': (7, 7, 7),
                    'gamma': True
                },
                txt='gs_licoo2.txt')

    a.set_calculator(calc)
    a.get_potential_energy()
    parprint('{} Saving calculation \n'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    calc.write('licoo2.gpw', mode='all')

    parprint('{} Band structure calc \n'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    calc = GPAW('licoo2.gpw',
                nbands=-20,
                fixdensity=True,
                symmetry='off',
                kpts={
                    'path': 'LGX',
                    'npoints': 120
                },
                convergence={'bands': -20},
                txt='licoo2_bands.txt')
    calc.get_potential_energy()
    calc.write('licoo2_bands.gpw', mode='all')

    bs = calc.band_structure()
    bs.plot(filename='bandstructure.png', show=False, emax=10.0, emin=-5.0)

    parprint('{} Saving band structure \n'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    with open('band_energies.pickle', 'wb') as handle:
        pickle.dump(calc.band_structure().todict(),
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL)
    parprint("{} Done".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

if False:
    seed = 'licoo2'

    calc = GPAW(seed + '.gpw', txt=None)

    w90.write_input(calc,
                    orbitals_ai=[[], [0, 1, 2, 3]],
                    bands=range(4),
                    seed=seed,
                    num_iter=1000,
                    plot=True)
    w90.write_wavefunctions(calc, seed=seed)
    os.system('wannier90.x -pp ' + seed)

    w90.write_projections(calc, orbitals_ai=[[], [0, 1, 2, 3]], seed=seed)
    w90.write_eigenvalues(calc, seed=seed)
    w90.write_overlaps(calc, seed=seed)