from ase.io import read
from gpaw import GPAW, PW, FermiDirac
import os
licoo2 = read("small.cif")

if os.path.isfile('licoo2_gs.gpw') == False:
    print("restart file not found :(")
    calc = GPAW(mode=PW(600),
                xc='PBE',
                kpts=(8, 8, 8),
                # random guess (needed if many empty bands required)
                random=True,
                occupations=FermiDirac(0.01),
                nbands=26,
                txt='outputs/licoo2_gs.txt')
    licoo2.calc = calc
    licoo2.get_potential_energy()
    calc.write('licoo2_gs.gpw')

calc = GPAW('licoo2_gs.gpw',
            nbands=30,
            fixdensity=True,
            symmetry='off',
            kpts={'path': licoo2.cell.bandpath().path, 'npoints': 60},
            convergence={'bands': 20},
            txt='outputs/licoo2_bands.txt')
calc.get_potential_energy()
calc.write('licoo2_bands.gpw', mode='all')
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=False, emax=15.0,
        emin=-4.0)
