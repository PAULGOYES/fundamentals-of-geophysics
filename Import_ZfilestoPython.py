#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:31:48 2020

@author: Paul Goyes
This Script imports you Z_file from Stratagem HS4 used for MT sounding
E-mail: goyes.yesid@gmail.com
E-mail: ypgoype@uis.edu.co
Universidad Industrial de Santander
Geological Faculty
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
filepath = 'ZFILE.001'

data = []

headers = ['frequency',
        'ExHy_coherency',
        'ExHy_scalar_apparent_resistivity',
        'ExHy_scalar_phase',
        'EyHz_coherency',
        'EyHx_scalar_apparent_resistivity',
        'EyHx_scalar_phase',
        're_Zxx/sqrt(µo)',
        'im_Zxx/√(µo)',
        're_Zxy/√(µo)',
        'im_Zxy/√(µo)',
        're_Zyx/√(µo)',
        'im_Zyx/√(µo)',
        're_Zyy/√(µo)',
        'im_Zyy/√(µo)',
        ]



with open(filepath) as fp:
    while True:
        line = fp.readline()
        if not len(line):
            break

        fp.readline()
        line2 = fp.readline()
        fp.readline()

        combined = line.strip().split() + line2.strip().split()

        data.append(combined)
        

df = pd.DataFrame(data, columns=headers).astype('float')
array = np.array(data).astype(np.float)

# example of type
#print(type(df['frequency'][0]))

xx = df.frequency.values
yy = df.ExHy_scalar_apparent_resistivity.values
phase = df.ExHy_scalar_phase.values

plt.subplot(211)
plt.loglog(xx[1:40],yy[1:40],'b-o')
plt.xlabel('Frequency [Hz]')
plt.ylabel('ExHy_scalar_apparent_resistivity')
plt.grid(True)
plt.subplot(212)
plt.semilogx(xx[1:40],phase[1:40],'b-o')
plt.xlabel('Frequency [Hz]')
plt.ylabel('ExHy_scalar_phase')
plt.grid(True)
