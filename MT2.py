#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 21:44:04 2020

@author: dharma
"""

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from SimPEG import Mesh, Utils, Solver
from scipy.constants import mu_0, epsilon_0
import matplotlib


from MTforward import simulateMT
# help(simulateMT)

rho_half = 100.  # Resistivity of the halfspace in Ohm-m
sigma_half = 1./rho_half  # Conductivity is the inverse of resistivity

frequency = np.logspace(-3, 2, 25)  # frequencies at which to simulate the MT problem

##############
def skin_depth(sigma, f):
    """
    Depth at which the fields propagating through a homogeneous 
    medium have decayed by a factor of 1/e for a given 
    frequency, f and conductivity, sigma
    """
    return 500./np.sqrt(sigma * f)

skin_depth_min = skin_depth(sigma_half, frequency.max())
skin_depth_max = skin_depth(sigma_half, frequency.min())

print("The minimum skin depth is {:1.2f}m".format(skin_depth_min))
print("The maximum skin depth is {:1.2e}m".format(skin_depth_max))
#####################

cs = skin_depth_min / 4.
core_extent = 5000. 
domain_extent = 2 * skin_depth_max

print("The smallest cell size is {:1.2f}m".format(cs))
print("The core region of the mesh extends {:1.2e}m".format(core_extent))
print("The mesh should extend at least {:1.2e}m".format(domain_extent))

npad = 1  # start with 1 cell
padding_fact = 1.3  # the amount by which we will expand each cell of the padding

def padding_extent(npad):
    """
    given a number of padding cells, this computes how far 
    the padding extends
    """
    padding_widths = cs*padding_fact**(np.arange(npad) + 1)
    return padding_widths.sum()

# keep adding padding until we are beyond the desired extent
padding_z = padding_extent(npad)
while padding_z < domain_extent:
    npad+=1
    padding_z = padding_extent(npad)
    
print(
    "{:1.0f} padding cells extends {:1.2e}m > {:1.2e}m "
    "(2 skin depths)".format(
        npad, padding_extent(npad), domain_extent
    )
)
    
ncz = np.ceil(core_extent / cs)  # number of cells in the core domain
hz = [(cs, npad, -1.3), (cs, ncz)]  # define how to construct the cell widths
mesh = Mesh.TensorMesh([hz], x0='N')  # construct a 1D Tensor Mesh

print(
    "There are {:1.0f} cells in the mesh. The mest extends {:1.2e}m".format(
        ncz, mesh.hx.sum()
    )
)


# plot the mesh
fig, ax = plt.subplots(1,1, figsize=(8, 3))
mesh.plotGrid(centers=True, faces=True, ax=ax)
ax.legend(["centers", "faces"])
ax.grid(which="both", linewidth=0.5)
ax.invert_xaxis()  # so that the surface is on our left hand side
ax.set_xlabel('z (m)')

# Set up a model

rho_target = 10.  # resistivity in Ohm-m
depth = 2000.  # depth to the top of the target in m 
thickness = 1000.  # thickness of the target in m


# put the model on the mesh
sigma = 1./rho_half * np.ones(mesh.nC)

# find the indices of the layer
layer_inds = (
    (mesh.vectorCCx<=-depth) & 
    (mesh.vectorCCx>-(depth+thickness))
)
sigma[layer_inds] = 1./rho_target

# plot the model
fig, ax = plt.subplots(1, 1, figsize=(8, 3))

# trickery to plot from node to node rather than at cell centers
z = np.repeat(mesh.vectorNx[1:-1], 2, axis=0)
z = np.r_[mesh.vectorNx[0], z, mesh.vectorNx[-1]]
sigma_plt = np.repeat(sigma, 2, axis=0)

ax.semilogy(z, sigma_plt,"C0", lw=2)
ax.grid(which="both", linewidth=0.5)
ax.set_xlim([-5000., 0.])
ax.set_ylim([5e-3, 1])
ax.invert_xaxis() # plot the surface on the left

ax.set_xlabel("Elevation (m)", fontsize=14)
ax.set_ylabel("Conductivity (S/m)", fontsize=14)


from SimPEG.EM.Analytics import MT_LayeredEarth


# the analytic takes the frequencies, layer thicknesses and layer conductivities
sigma_layers = np.r_[
    1./rho_half, 
    1./rho_target, 
    1./rho_half]
h = np.r_[depth, thickness]  

app_res_ana, app_phase_ana = MT_LayeredEarth(
    frequency, h, sigma_layers, 'Res-Phase'
)
# numerically compute the response
app_res, app_phase = simulateMT(mesh, sigma, frequency)

def plot_with_analytic(frequency, app_res, app_phase, app_res_ana, app_phase_ana):
    # Plot and compare the results
    fig, ax = plt.subplots(2, 1, figsize=(8, 3*2))

    # apparent resistivity
    ax[0].loglog(frequency, app_res, label='Numeric')
    ax[0].loglog(frequency, app_res_ana, 'k.', label='Analytic')
    ax[0].set_ylabel("$\\rho_a \ (\Omega m)$", fontsize=14)
    ax[0].set_ylim([2e1, 3e2])

    # phase
    ax[1].semilogx(frequency, app_phase, label='Numeric')
    ax[1].semilogx(frequency, app_phase_ana, 'k.', label='Analytic')
    ax[1].set_ylabel("$\phi \ (^{\circ})$", fontsize=13)
    ax[1].set_ylim([0., 90.])

    for a in ax:
        a.grid(True, which='both', linewidth=0.3)
        a.set_xlim(frequency.max(), frequency.min())
        a.set_xlabel("Frequency (Hz)", fontsize=14)
        a.legend(fontsize=11)

    plt.tight_layout()
    plt.show()
    
plot_with_analytic(frequency, app_res, app_phase, app_res_ana, app_phase_ana)

