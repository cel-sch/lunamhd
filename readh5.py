# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:59:04 2023

@author: celin
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from VenusMHDpy.Bsplines import Bspline
from VenusMHDpy.Grid import GRID

grid = GRID() 

def ploth5(file, varnr=0):
    """
    Reads h5py file for given eigenvalue run. Plots the eigenfunctions.
    """
    with h5py.File(file, 'r') as f:
        ms = f['grid']['m'][()] # array of poloidal mode numbers
        nus = f['grid']['nu'][()] # grid.nu is determined by which model is being used, for IdealMHDFlow-Euler is [3,2,2,2,3,3,2,2] because 8 variables
        N = f['grid']['N'][()] # number of grid points?
        #S = f['grid']['S'][()]
        
        ### Update grid with new parameters
        grid.N = N
        grid.nu = nus
        grid.S = f['grid']['S'][()]
        
        grid.Mmin = min(f['grid']['m'][()])
        grid.Mmax = max(f['grid']['m'][()])
        
        grid.knots = f['grid']['knots'][()]
        grid.sk = f['grid']['sk'][()]
        
        # not all grid parameters are specified here, bspline uses: knots, sk, N, S
        
        ### Build Bspline arrays
        r = np.linspace(0.,1.,10000)
        BspCalc = []
        
        nu = nus[varnr]
        fig, ax1 = plt.subplots()
        
        vec_allms = f['variables'][f'var{varnr}'] # all poloidal mode numbers for variable varnr (in h5py vars go from 0 to 7)

        #Create the Bspline arrays
        #-------------------------------------------------------------------------------------------------
        if nu not in BspCalc:
            BspCalc.append(nu)	
            vars()['Bsp'+str(nu)] = np.zeros(shape=(r.size, N+1+nu))
            
            for j in range(N):
                l = np.ones(len(r), dtype=int)*j-nu
                vars()['Bsp'+str(nu)][:,j] = Bspline(r,l,nu,grid,der=0)
                
        mixB = max(nus)-nu # 8 - order of variable (2 or 3 generally)
        # labels
        for j in range(f['grid']['Mtot'][()]): # selects vector for every poloidal mode number
            vec = vec_allms[j]
            mode = np.sum(vars()['Bsp'+str(nu)]*vec, axis=1) # performs a bspline on the vector?
            ax1.plot(r, mode, label=f'm={ms[j]}')
            
        ax1.legend(loc='best')
        ax1.set_title(f'var {varnr + 1}, file={file}')     
         
    