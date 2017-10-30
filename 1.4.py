from __future__ import print_function
import numpy as np
from ase import Atoms
from ase.units import eV, Ang, GPa, kJ
import Morse
import matplotlib.pyplot as plt
from ase.build import bulk

cu = bulk("Cu", "fcc", a=3.6, cubic=True)
calc = Morse.MorsePotential()
cu.set_calculator(calc)
#apply shear strain
cell = cu.get_cell()

x_strain = 0.1

n=10

stress = [10]
strain = np.linspace (-0.1, 0.1, n)

def poi(strain, stress, x_strain):
    cell[0,0] = 3.6*(1+x_strain)
    k=0
    if stress[0] > (0+1e-6) or stress[0] < (0-1e-6):
        stress=[]
        strain_calc = []
        for i in strain:
            cell[1,1]=cell[2,2]=3.6*(1+i)
            cu.set_cell(cell, scale_atoms=True)
            stress.append(cu.get_stress(voigt=False)[1,1])
            strain_calc.append(i)
            k += 1
            if k>1:
                if stress[-2]*stress[-1] < 0:
                    break
        strain = np.linspace (strain_calc[-2], strain_calc[-1],n)
        return poi(strain, stress, x_strain)
    else:
        min_strain = strain[0]
        return min_strain

p=50
x_strain = np.linspace (-0.1, 0.1, p)
poisson = []

for x in x_strain:
    #print(p)
    p -= 1
    min_strain = poi(strain, stress, x)
    poisson.append (abs(min_strain/x))
    
plt.plot(x_strain, poisson)
plt.show();
    
#min_strain = poi(strain, stress, x_strain)
#print('Minimum strain:', min_strain)
#print("Poisson's ratio:", )