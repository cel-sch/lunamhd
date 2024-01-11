# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 17:28:43 2023

@author: celin
"""

import matplotlib.pyplot as plt
import numpy as np

### VARYING ROTATION? ###
vrot = True

### READ LUNA DATA ###
lunarun = '0xc696d'
lunafile = np.load(f'Output/KH/{lunarun}/{lunarun}.npz', allow_pickle = True)

params = lunafile['params'].item()
profparams = lunafile['profparams'].item()
eps_a = params['r0']/params['R0']
beta = profparams['beta0']/eps_a**2

lunaevs = lunafile['eigenvals'] # changes normalisation
lunax = lunafile['scanvals']
if vrot:
    lunaomega = [i*np.sqrt(beta) for i in lunax]
    lunax = lunaomega

### READ AE DATA ###
AErun = '25-10-23_19;08_rhoP_Omega'
AEfile = np.load(f'AEs/{AErun}.npz', allow_pickle = True)

AEevs = AEfile['eigenvals'] # changes normalisation
#AEasyevs = [i.real*eps_a for i in AEfile['asy_eigenvals']]
AEx = AEfile['scanvals'] # for varying rotation: this is Omega

### GROWTH RATES ###
lunagams = [i.real for i in lunaevs]
AEgams = [i.imag*eps_a for i in AEevs]


### PLOT THINGS ###
if vrot:
    xlabel = '$Ω/(ε_aω_A)$'

fig, ax = plt.subplots()
ax.set_xlabel(f'{xlabel}',fontsize=12)
ax.set_ylabel('$γ/(ε_aω_A)$',fontsize=12)
ax.grid()

ax.plot(AEx[15:], AEgams[15:], '-x', label='$\hat{Ω}$ step: step model')
ax.plot(lunax[:-1], lunagams[:-1], '-x', label='$\hat{Ω}$ step: VENUS-MHD')
plt.text(0.5,0.076,'$(1,1)$',fontsize=16)
plt.text(0.5,0.068,'$\hat{β}=1.8$',fontsize=16)
plt.text(0.5,0.06,'$\hat{Δq}=1$',fontsize=16)
#plt.text(1.1,0.02,'$Δr_0=0.25$',fontsize=16)

ax.legend()

fig.set_size_inches(7, 5.5)