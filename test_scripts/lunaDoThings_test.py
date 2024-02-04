# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:11:19 2023

@author: celin
"""

import lunaSolverKH_test as lunsol
import readh5 as rh5
import matplotlib.pyplot as plt
import numpy as np
import os

### BEGIN ###
sol = lunsol.Solver()

### INITIALISE NEW VMEC EQUILIBRIUM? ###
sol.initialisers['RunVMEC'] = True
if sol.initialisers['RunVMEC'] == False:
    sol.VMEClabel = 'oldVMECrun' # if using an existing VMEC run, name of old run

os.system(f'mkdir -p Output/{sol.VMEClabel}') 
sol.out_filepath = f'Output/{sol.VMEClabel}'
sol.in_filename = 'rhoP_test_KH.in'

runVENUS = True
if runVENUS == False:
    out_filename = '0x4a2257.npz'
    out_file = f'{sol.out_filepath}/{out_filename}'
else:
    out_file = None  

### SET INITIALISERS ###
sol.initialisers['EVg_type'] = 'polynom_EV' # pick eigenvalue guess method
sol.initialisers['ToPlot'] = True # plot the VMEC profiles, EV and EF
# NOTE: setting ToPlot to True will interrupt the VENUS scan?

### LOAD DATA ###
data = sol.getData(dataFile = out_file, runSol = runVENUS)

### PLOT DATA ###
plotEigs = True
plotEigGuess = False
plotEigFuncs = True # plots first and last eigenfunction from scan for sanity checking

if plotEigs:
    ws = data['eigenvals']
    gams = [i.real for i in ws]
    if data['scanparams'] == 'mach':
        #x = data['v0vas']
        #xlabel = 'v0/va'
        profparams = data['profparams'].item()
        params = data['params'].item()
        
        eps_a = params['r0']/params['R0']
        betahat = profparams['beta0']/eps_a**2 # normalised beta
        mach = data['scanvals']
        Omega = [i*np.sqrt(betahat) for i in mach]  # normalised omega
        
        x = Omega
        xlabel = 'Omega'
    else:
        x = data['scanvals']
        xlabel = data['scanparams']
    
    fig, ax = plt.subplots()
    ax.plot(x[:], gams[:], '-x',label='VENUS-MHD $\hat{γ}$')
    if plotEigGuess:
        wsguess = data['eigenguesses']
        gamguess = [i.real for i in wsguess]
        ax.plot(x, gamguess, '-x', label='EV guess')
    ax.set_ylabel('gamma')
    ax.set_xlabel(f'{xlabel}')
    ax.set_title(f'{sol.VMEClabel}') # not perfect
    if runVENUS == False:
        ax.set_title(f'{out_filename}'.replace('.npz',''))
    plt.grid()
    if plotEigGuess:
        plt.legend()
        
if plotEigFuncs: 
    Nscan = len(data['scanvals'])
    label = sol.VMEClabel
    
    if Nscan == 0:
        filename = f'{label}_0.hdf5'
        file = f'{sol.out_filepath}/{filename}'
        rh5.ploth5(file, allVars=True)
        
    else:
        filename0 = f'{label}_0.hdf5'
        file0 = f'{sol.out_filepath}/{filename0}'
        
        filenameN = f'{label}_{Nscan-1}.hdf5'
        fileN = f'{sol.out_filepath}/{filenameN}'
        
        rh5.ploth5(file0)
        rh5.ploth5(fileN)
    
plt.show()
    

### GET SINGLE RUN PLOTS ###
# Get plots with profiles and eigenfunctions etc
# Need to find a way to do this without re-running the eigenfunction guess
# Also this loads default parameters for IdealMHDFlow-Euler model?
plotSingle = False
if plotSingle:
    sol.VMEClabel = '0x73ba54'
    sol.initialisers['ToPlot'] = True
    sol._runVENUS(labelnr=0)
    

