# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 11:46:00 2023

@author: celin
"""
import matplotlib.pyplot as plt
import numpy as np

### VARYING ROTATION? ###
vrot = True

### READ LUNA DATA ###
#lunarun = '0x1197a5'
lunarun = '0x70960f'
lunafile = np.load(f'Output/KH/{lunarun}/{lunarun}.npz', allow_pickle = True)

lunarun2 = '0xb4bb63'
lunafile2 = np.load(f'Output/KH/{lunarun2}/{lunarun2}.npz', allow_pickle = True)

params = lunafile['params'].item()
eps_a = params['r0']/params['R0']

profparams = lunafile['profparams'].item()
beta = profparams['beta0']/eps_a**2

params2 = lunafile2['params'].item()
profparams2 = lunafile2['profparams'].item()
beta2 = profparams2['beta0']/eps_a**2

###

lunaevs = lunafile['eigenvals'] # changes normalisation
lunax = lunafile['scanvals']

lunaevs2 = lunafile2['eigenvals'] # changes normalisation
lunax2 = lunafile2['scanvals']
if vrot:
    lunaomega = [i*np.sqrt(beta) for i in lunax]
    lunaomega2 = [i*np.sqrt(beta2) for i in lunax2]
    lunax = lunaomega
    lunax2 = lunaomega2

### READ AE DATA ###
AErun = '24-10-23_17;42_rhoP_Omega'
AEfile = np.load(f'AEs/{AErun}.npz', allow_pickle = True)

AEevs = AEfile['eigenvals'] # changes normalisation
#AEasyevs = [i.real*eps_a for i in AEfile['asy_eigenvals']]
AEx = AEfile['scanvals'] # for varying rotation: this is Omega

### GROWTH RATES ###
lunagams = [i.real for i in lunaevs]
lunagams2 = [i.real for i in lunaevs2]
AEgams = [i.imag*eps_a for i in AEevs]


### PLOT THINGS ###
if vrot:
    xlabel = '$Ω/(ε_aω_A)$'

fig, ax = plt.subplots()
ax.set_xlabel(f'{xlabel}',fontsize=12)
ax.set_ylabel('$γ/(ε_aω_A)$',fontsize=12)
ax.grid()

ax.plot(AEx[:], AEgams[:], '-x', label='$ρ, p$ step: step model')
ax.plot(lunax[:], lunagams[:], '-x', label='$ρ, p$ step: VENUS-MHD')
#ax.plot(lunax2[:], lunagams2[:], '-x', label='$ρ, p$ step: VENUS-MHD')
ax.legend()

plt.text(1,0.9,'$(2,1)$',fontsize=16)
plt.text(1,0.8,'$\hat{β}=0.5$',fontsize=16)
plt.text(1,0.7,'$\hat{Δq}=1$',fontsize=16)
#plt.text(1.45,0.075,'$Δr_0=0.4$',fontsize=16)
fig.set_size_inches(7, 5.5)
