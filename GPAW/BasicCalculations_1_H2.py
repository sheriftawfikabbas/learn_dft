from ase import Atoms
from gpaw import GPAW

# Setting up the calculation

# Place one H atom at point (5,5,5) in a crystal of size 10 x 10 x 10:

atoms = Atoms('H', positions=[(5, 5, 5)], cell=(10, 10, 10))

# Place one H2 molecule, where one H is at (5-0.35,5,5) and the other is at (5+0.35,5,5) in a crystal of size 10 x 10 x 10. Here, the list of positions of atoms go into positions. Note that: here, positions becomes a list of tuples

atoms = Atoms('H2', positions=[(5-0.35, 5, 5),
                               (5+0.35, 5, 5)], cell=(10, 10, 10))

# and it can also be a list of lists

atoms = Atoms('H2', positions=[[5-0.35, 5, 5],
                               [5+0.35, 5, 5]], cell=(10, 10, 10))

# Let' use functions to get properties of the system

atoms.get_positions()

atoms.get_chemical_symbols()

atoms.get_distance(0, 1)

atoms.get_initial_charges()

atoms.get_initial_magnetic_moments()

atoms.get_masses()

atoms.get_number_of_atoms()

atoms.get_pbc()

atoms.get_volume()

atoms.get_cell_lengths_and_angles()

# Now let's use variables for more compact notation

a = 10
d = 0.74

atoms = Atoms('H2', positions=[[a/2-d/2, a/2, a/2],
                               [a/2+d/2, a/2, a/2]], cell=(a, a, a))

# Next, we need to create a calculator object

# Calculation of the atomization energy
# This is the energy of splitting all atoms in the molecule into isolated atoms.
# So it quantifies the stength of the bond(s) in the molecule.

# First let's calculate the total energy of a single H atom

H = Atoms('H', positions=[(a/2, a/2, a/2)], cell=(a, a, a))

# Then create the calculator

# gpaw calculator: create the GPAW object and then set this object as the calculator for the atoms object
calc = GPAW()
H.set_calculator(calc)
# Then, calculate the total energy. That's the first DFT calculation we will do.

e = H.get_potential_energy()
print(e)

# Second, we calculate the total energy of H2

H2 = Atoms('H2', positions=[[a/2-d/2, a/2, a/2],
                            [a/2+d/2, a/2, a/2]], cell=(a, a, a))
H2.set_calculator(calc)
eTotal = H2.get_potential_energy()
print(eTotal)

# The atomization energy can be calculated as: eTotal - 2 x e

print(eTotal-2*e)

# But that's pretty large, right? H2 can be broken by ~ -4.5 eV, not -6.6 eV!

# How do we fix this?

# Check this out: the spin of the single H atom is zero, but this doesn't make sense!
H.get_magnetic_moments()

# To fix this, we need to apply Hund's rule on the H atom

# Note: Hund rule setting is only valid for single atoms

calc_HundRule = GPAW(hund=True)
H.set_calculator(calc_HundRule)
e = H.get_potential_energy()
print(e)

# The value of -0.876647212841809 eV is much different for sure.

print(eTotal-2*e)

# Now the value is much closer to the correct value: it is ~ -4.9 eV.

# There is still room for improvement: the exchange-correlation. Let's use PBE.

calc_HundRule_PBE = GPAW(hund=True, xc='PBE')
H.set_calculator(calc_HundRule_PBE)
e = H.get_potential_energy()
calc_PBE = GPAW(xc='PBE')
H2.set_calculator(calc_PBE)
eTotal = H2.get_potential_energy()
print(eTotal-2*e)

# Now, -4.488165097141112 is much closer to -4.5 eV.