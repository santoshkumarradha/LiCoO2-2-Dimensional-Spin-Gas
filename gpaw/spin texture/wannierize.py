from gpaw import restart
from ase.dft import Wannier
import gpaw.wannier90 as w90
from gpaw import GPAW
import os

seed = 'licoo2_gs'

print("restarting and saving...")
calc = GPAW(seed + '.gpw', txt=None)
print("Done\n")

print("writing wannier input ... \n\n")
w90.write_input(calc, orbitals_ai=[[0, 1, 2, 3], [4, 5, 6, 7, 8], [], []],
                bands=[8, 9, 10, 11, 12, 13, 14, 15, 16],
                seed=seed,
                num_iter=1000,
                plot=True)
print("Done\n")


print("writing wannier wav ... \n\n")
w90.write_wavefunctions(calc, seed=seed)
print("Done\n")
sucess = 0
if sucess:
    print("running wannier 90...\n")
    os.system('wannier90.x -pp ' + seed)
    print("done\n")

    w90.write_projections(calc, orbitals_ai=[[0, 1, 2, 3], [
        4, 5, 6, 7, 8], [], []], seed=seed)
    w90.write_eigenvalues(calc, seed=seed)
    w90.write_overlaps(calc, seed=seed)
    os.system('wannier90.x ' + seed)
