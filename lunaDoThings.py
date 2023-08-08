# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:11:19 2023

@author: celin
"""

import lunaSolverKH as lunsol
import readh5 as rh5
import matplotlib.pyplot as plt
import numpy as np
import os

#plt.ion()

###
sol = lunsol.Solver()

### CALCULATE NEW EIGENVALUES OR NOT? ###
sol.initialisers['RunVMEC'] = False
if sol.initialisers['RunVMEC'] == False:
    sol.VMEClabel = '0xe5eb62'

os.system(f'mkdir -p Output/KH/{sol.VMEClabel}') 
#os.system(f'mkdir -p Output/KH')

sol.out_filepath = f'Output/KH/{sol.VMEClabel}'
#sol.out_filepath = f'Output/KH'
sol.in_filename = 'input_KH.in'

runVENUS = False # if VENUS is run without VMEC, I think the same run name is used so old .npz files will be replaced unless VMEClabel is changed?
if runVENUS == False:
    out_filename = '0xe5eb62.npz'
else:
    out_filename = 'blah.npz'    
out_file = f'{sol.out_filepath}/{out_filename}'


sol.initialisers['ToPlot'] = False

### LOAD DATA ###
data = sol.getData(dataFile = out_file, runSol = runVENUS)

### PLOT DATA ###
fig5 = False
AEs = True
plotEigs = True
plotEigGuess = False
plotEigFuncs = False

if plotEigs:
    ws = data['eigenvals']
    gams = [i.real for i in ws]
    if data['scanparams'] == 'mach':
        #x = data['v0vas']
        #xlabel = 'v0/va'
        profparams = data['profparams'].item()
        params = data['params'].item()
        
        eps_a = params['r0']/params['R0']
        beta = profparams['beta0']/eps_a**2
        mach = data['scanvals']
        Omega = [i*np.sqrt(beta) for i in mach]
        
        x = Omega
        xlabel = 'Omega'
        if fig5: ### Figure 5 ###
            mach5, gam5 = np.loadtxt('AEs/gam_fig5.txt', dtype=(np.cfloat)).transpose()
            gam5 = [i.real*eps_a for i in gam5] # eps_a normalisation removed to match VENUS normalisation
            Omega5 = [i.real*np.sqrt(beta) for i in mach5]
        if AEs:
            aeName = '0xbd1b8b'
            aeFile = np.load(f'AEs/{aeName}.npz', allow_pickle = True)
            
            aegams = [i.imag*eps_a for i in aeFile['eigenvals']]
            a_aegams = [i.imag*eps_a for i in aeFile['asy_eigenvals']]
            aex = aeFile['scanvals']
    else:
        x = data['scanvals']
        xlabel = data['scanparams']
    
    fig, ax = plt.subplots()
    ax.plot(x[:], gams[:], '-x',label='VENUS-MHD $\hat{γ}$')
    if fig5:
        ax.plot(Omega5[:13], gam5[:13], '-.', label='2013 PPCF $\hat{γ}$')
    if AEs:
        ax.plot(aex[23:], aegams[23:], '-.', label='step model $\hat{γ}$')
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
    if plotEigGuess or fig5 or AEs:
        plt.legend()
        
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
    sol.VMEClabel = '0xf3cff5'
    sol.initialisers['ToPlot'] = True
    sol._runVENUS(labelnr=0)
    

