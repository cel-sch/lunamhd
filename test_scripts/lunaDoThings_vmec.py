# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:11:19 2023

@author: celin
"""
# Script for producing VMEC input files

import lunaSolverKH as lunsol
import readh5 as rh5
import matplotlib.pyplot as plt
import numpy as np
import os

### BEGIN ###
sol = lunsol.Solver()

### CALCULATE NEW EIGENVALUES OR NOT? ###
sol.initialisers['RunVMEC'] = True

os.system(f'mkdir -p Output/vmec_testing/{sol.VMEClabel}') 
#os.system(f'mkdir -p Output/KH')

sol.out_filepath = f'Output/vmec_testing/{sol.VMEClabel}'
sol.in_filename = 'vmec_test.in'

runVENUS = True # if VENUS is run without VMEC, I think the same run name is used so old .npz files will be replaced unless VMEClabel is changed?
if runVENUS == False:
    out_filename = '0xe5eb62.npz'
else:
    out_filename = 'blah.npz'    
out_file = f'{sol.out_filepath}/{out_filename}'


### SET INITIALISERS ###
sol.initialisers['ToPlot'] = True
sol.initialisers['EVg_type'] = 'polynom_EV'

### LOAD DATA ###
data = sol.getData(dataFile = out_file, runSol = runVENUS)

### PLOT DATA ###
plotEigs = True
plotEigFuncs = True

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
    ax.set_ylabel('gamma')
    ax.set_xlabel(f'{xlabel}')
    ax.set_title(f'{sol.VMEClabel}') # not perfect
    if runVENUS == False:
        ax.set_title(f'{out_filename}'.replace('.npz',''))
    plt.grid()
        
if plotEigFuncs: # what is input name if VMEC not run but VENUS is?
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
    

