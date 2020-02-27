#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Feb 27 11:31:48 2020

@author: Modified by Paul Goyes
E-mail: goyes.yesid@gmail.com
E-mail: ypgoype@uis.edu.co
Universidad Industrial de Santander

THIS SCRIPT CREATES THE FIG 1 OF
http://dx.doi.org/10.1190/tle36080696.1
"""

'''
====================================
1D MAGNETOTELLURIC MODELLING PROGRAM
====================================
FIRST UPDATED 17TH DECEMBER 2013
ORIGINAL VERSION DEVELOPED BY ANDREW PETHICK   
WWW.DIGITIALEARTHLAB.COM      
'===================================='''

import math
import cmath
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

start = time.clock();

# Please check the Magnetic Permeability (H/m)
# at: https://docs.scipy.org/doc/scipy-0.14.0/reference/constants.html#module-scipy.constants
mu = constants.mu_0; #Magnetic Permeability (H/m)

# Parameters of our resistivity-model
resistivities = [1.0e+02,1.0e+01,1.0e+02]; # in Ohm-meter
thicknesses = [2000, 1000]; # In meters
frequencies=np.logspace(-3, 2, 25) # In logarithmic equidistant range. (Hz)

# forward modelling operator F[m] = d
n = len(resistivities);
apparentResistivity=[]
phase=[]

for frequency in frequencies:   
    w =  2*math.pi*frequency;       
    impedances = list(range(n));
    #compute basement impedance
    impedances[n-1] = cmath.sqrt(w*mu*resistivities[n-1]*1j);

    for j in range(n-2,-1,-1):
        resistivity = resistivities[j];
        thickness = thicknesses[j];
  
        # 3. Compute apparent resistivity from top layer impedance
        # Step 2. Iterate from bottom layer to top(not the basement) 
        # Step 2.1 Calculate the intrinsic impedance of current layer
        dj = cmath.sqrt((w * mu * (1.0/resistivity))*1j);
        wj = dj * resistivity;
        # Step 2.2 Calculate Exponential factor from intrinsic impedance
        ej = cmath.exp(-2*thickness*dj);                     
    
        # Step 2.3 Calculate reflection coeficient using current layer
        #          intrinsic impedance and the below layer impedance
        belowImpedance = impedances[j + 1];
        rj = (wj - belowImpedance)/(wj + belowImpedance);
        re = rj*ej; 
        Zj = wj * ((1 - re)/(1 + re));
        impedances[j] = Zj;    

    # Step 3. Compute apparent resistivity from top layer impedance
    Z = impedances[0];
    absZ = abs(Z);
    apparentResistivity.append((absZ * absZ)/(mu * w))
    phase.append(math.atan2(Z.imag, Z.real))
    
print('');
print('time taken = ', time.clock() - start, 's');

# Plot the results in 2 axes. Phase (degree) and apparent resistivity (ohm meter)

fig, ax = plt.subplots(2, 1, figsize=(6, 4))

ax[0].loglog(frequencies,apparentResistivity,'r-o'), 
ax[0].set_ylabel("$\\rho_a \ (\Omega m)$", fontsize=12)

ax[1].semilogx(frequencies,np.array(phase)*180/np.pi,'b-o'),
ax[1].set_ylabel("$\phi \ (^{\circ})$", fontsize=12)
ax[1].set_ylim([0., 90.])

for a in ax:
    a.grid(True, which='both', linewidth=0.2)
    a.set_xlim(frequencies.max(), frequencies.min())
    a.set_xlabel("Frequency (Hz)", fontsize=10)

plt.tight_layout()
plt.show()
