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

grid = GRID() # initialises a default grid

#out_filepath = 'Output/KH/0x24b7c7'
#filename = '0x24b7c7_5.hdf5'

#file = f'{out_filepath}/{filename}'

def ploth5(file, allVars=False):
    """
    Reads h5py file for given eigenvalue run. Plots the eigenfunctions.
    
    allVars - If True, plots every eigenfunction (vars 1-8). Otherwise, only
        var 1 is plotted.
    """
    with h5py.File(file, 'r') as f:
        ms = f['grid']['m'][()] # array of poloidal mode numbers
        nus = f['grid']['nu'][()]
        N = f['grid']['N'][()]
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
        
        if allVars:
            for i,nu in enumerate(nus): # grid.nu is determined by which model is being used, for IdealMHDFlow-Euler is [3,2,2,2,3,3,2,2] because 8 variables		
            # counting from 1 to 8
                fig, vars()['ax'+f'{i}'] = plt.subplots()
                
                var = list(f['variables'])[i]
                vec_allms = f['variables'][f'{var}'] # all poloidal mode numbers for a given variable

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
                    vars()['ax'+f'{i}'].plot(r, mode, label=f'm={ms[j]}')
                    
                vars()['ax'+f'{i}'].legend()
                vars()['ax'+f'{i}'].set_title(f'var {i+1}, file={file}')
        else:
            nu = nus[0]
            fig, ax1 = plt.subplots()
            
            vec_allms = f['variables']['var0'] # all poloidal mode numbers for var 1 (in h5py vars go from 0 to 7)

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
                
            ax1.legend()
            ax1.set_title(f'var 1, file={file}')     
        
#ploth5(file)        
#plt.show()