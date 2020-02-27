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

start = time.clock();

mu = 4*math.pi*1E-7; #Magnetic Permeability (H/m)
#resistivities = [300, 2500, 0.8, 3000, 2500];
#thicknesses = [200, 400, 40, 500];
resistivities = [1.0e+02,1.0e+01,1.0e+02];
thicknesses = [2000, 1000];


#frequencies = [0.0001,0.005,0.01,0.05,0.1,0.5,1,5,10,50,100,500,10000];
frequencies=np.logspace(-3, 2, 25) 
#print('freq\tares\t\t\tphase');

n = len(resistivities);
# crea una lista vacia de resistividades aparentes
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
        #Step 2. Iterate from bottom layer to top(not the basement) 
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
#    print(frequency, '\t', apparentResistivity, '\t', phase);
    
print('');
print('time taken = ', time.clock() - start, 's');


import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 1, figsize=(8, 3*2))

ax[0].loglog(frequencies,apparentResistivity,'r-o'), 
ax[0].set_ylabel("$\\rho_a \ (\Omega m)$", fontsize=14)

ax[1].semilogx(frequencies,np.array(phase)*180/np.pi,'b-o'), plt.grid()
ax[1].set_ylabel("$\phi \ (^{\circ})$", fontsize=13)
ax[1].set_ylim([0., 90.])

for a in ax:
    a.grid(True, which='both', linewidth=0.3)
    a.set_xlim(frequencies.max(), frequencies.min())
    a.set_xlabel("Frequency (Hz)", fontsize=14)

plt.tight_layout()
