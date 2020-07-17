from ase.io import read
from gpaw import GPAW, PW, FermiDirac
import os
import pymatgen as p
from pymatgen.io import ase
import re
import numpy as np
from ase.parallel import parprint
from pymatgen.io.cif import CifParser

def get_ase(n=3,m=1):
    parser = CifParser("small.cif")
    lco = parser.get_structures()[0]
    lco_sc=lco.copy()
    n=3
    sc=np.array([
        [n,0,0],
        [0,n,0],
        [0,0,n]
    ])
    lco_sc.make_supercell(sc)
    li=np.array([i if j.species_string=="Li" else -1 for i,j in enumerate(lco_sc)])
    li=li[li>=0]
    index = np.random.choice(li.shape[0], m, replace=False)  
    for i in index:
        lco_sc.replace(li[i],"H")
    percentage_li=lco_sc.composition.get_atomic_fraction("Li")/lco_sc.composition.get_atomic_fraction("Co")
    # lco_sc.to("cif","lco_sc.cif")
    return ase.AseAtomsAdaptor().get_atoms(lco_sc),percentage_li

licoo2,percentage=get_ase()

parprint("percentage of defect concentration = {:.3f}".format(percentage))

if os.path.isfile('licoo2_gs.gpw') == False:
    print("restart file not found :(")
    calc = GPAW(mode=PW(600),
                xc='PBE',
                symmetry={'point_group': False},
                kpts=(3, 3, 3),
                # random guess (needed if many empty bands required)
                random=True,
                occupations=FermiDirac(0.01),
                nbands=-30,
                txt='outputs/licoo2_gs.txt')
    licoo2.calc = calc
    licoo2.get_potential_energy()
    calc.write('licoo2_gs.gpw')

calc = GPAW('licoo2_gs.gpw',
            nbands=-30,
            fixdensity=True,
            symmetry='off',
            kpts={'path': licoo2.cell.bandpath().path, 'npoints': 60},
            convergence={'bands': -30},
            txt='outputs/licoo2_bands.txt')
calc.get_potential_energy()
calc.write('licoo2_bands.gpw', mode='all')
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=False, emax=15.0,
        emin=-5.0)

