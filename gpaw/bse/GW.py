from gpaw.response.g0w0 import G0W0
from ase.build import mx2
from gpaw import GPAW, PW, FermiDirac


restart = 'licoo2_gs.gpw'

calc = GPAW(restart)
calc.diagonalize_full_hamiltonian()
calc.write('licoo2_fulldiag.gpw', 'all')
for ecut in [300]:
    gw = G0W0(calc='licoo2_fulldiag.gpw',
              bands=(10, 11),
              ecut=ecut,
              #   truncation='2D',
              nblocksmax=True,
              #   q0_correction=True,
              filename='licoo2_g0w0_{}'.format(ecut),
              savepckl=True)

    gw.calculate()
