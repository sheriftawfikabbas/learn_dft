from ase import Atoms
from gpaw import GPAW
from ase.optimize import QuasiNewton

# Setting up the calculation

# Place one CO2 molecule in the center of the cell

a = 10
d = 1.17103

CO2 = Atoms('CO2', positions=[[a/2,a/2,a/2],[a/2-d,a/2,a/2],[a/2+d,a/2,a/2]], cell=(a,a,a))

# Next, we need to create a calculator object

# The optimization calculation
# First we create the GPAW object

calc_PBE = GPAW(xc='PBE')
CO2.set_calculator(calc_PBE)

# By the way: by default, convergence of SCF is such that difference in energy is < 0.0005

# We can change this by adding the parameter: convergence={'energy': 0.0001}

# Next, we use a relaxation procedure
# Here, the atoms will be moved small steps in search for the optimal structure until the force on each atom is less that 0.05 eV/Ang

relax = QuasiNewton(CO2)
relax.run(fmax=0.05)
CO2.get_distance(0,1)

# We started from the correct structure, that's why the optimization finished in only 1 step.
# Change d to 1 and see what happens. The optimization will finish in a number of iterations.
