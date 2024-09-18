from ase import Atoms
from gpaw import GPAW
from ase.optimize import QuasiNewton
from ase.build import fcc111
from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms

# Setting up the calculation


slab = fcc111('Pt', size=(1, 1, 3))
slab.center(vacuum=10.0, axis=2)
slab.write('Pt_1x1x3.cif')


# Step 1: Calculate the total energy of the slab
# -----------------------------------------------

# It's s slab, so we need to "freeze" the atomic position of the the bottom of the slab.
# This is to simulate the sitation of a "thick" slab.
# First, let's find out what are the indices of the atoms:

slab.get_positions()

# So indeces 0 and 1 are the ones we want to fix

c = FixAtoms(indices=[0, 1])
slab.set_constraint(c)

# Like we did in the previous sheet: setup an optimization calculation using PBE.
# But now we are going to add spin polarization and kpoints.


calc_PBE = GPAW(xc='PBE', kpts=(2, 2, 1), spinpol=True)
slab.set_calculator(calc_PBE)
relax = QuasiNewton(slab)
relax.run(fmax=0.05)
Pt_slab_energy = slab.get_total_energy()

print("Slab energy: ", Pt_slab_energy)
# It is -20.505906792376535

slab.write('Pt_1x1x3_Optimized.cif')

# Step 2: Calculate the total energy of the slab + H
# -----------------------------------------------

add_adsorbate(slab, 'H', 1.7, 'ontop')
slab.write('Pt_1x1x3_Optimized_H.cif')

slab.set_calculator(calc_PBE)
relax = QuasiNewton(slab)
relax.run(fmax=0.05)
Pt_slab_H_energy = slab.get_total_energy()
slab.write('Pt_1x1x3_H_Optimized.cif')

print("Slab energy with H adsorbed: ", Pt_slab_H_energy)

# It is -24.228088785650773

# Recall from BasicCalculations_1_H2.py
H2_energy = -6.725007250781371

# The adsorption energy is:
E_ad = Pt_slab_H_energy - Pt_slab_energy - H2_energy/2

print("Adsorption energy:", E_ad)
# -0.39759964909709344

# To get closer to the ~ -0.5 eV result, we need more k-points