# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:11:19 2023

@author: celin
"""

import eigScannerIK as lunsol
import matplotlib.pyplot as plt
import numpy as np

#plt.ion()

sol = lunsol.Solver()

sol.out_filepath = 'ScanOutput/IK'
sol.in_filename = 'eigs_IK.in'

### CALCULATE NEW EIGENVALUES OR NOT? ###
sol.initialisers['RunVMEC'] = True
if sol.initialisers['RunVMEC'] == False:
    sol.VMEClabel = '0x7b3f1f' #'0xd9e513' for KH

runVENUS = True
if runVENUS == False:
    out_filename = '0xf3cff5.npz'
else:
    out_filename = 'blah.npz'    
out_file = f'{sol.out_filepath}/{out_filename}'


### PLOT DATA? ###
sol.initialisers['ToPlot'] = True


### LOAD DATA ###
data = sol.getData(dataFile = out_file, runSol = runVENUS)

