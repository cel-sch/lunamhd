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

def ploth5(file, varnr = 0):
    """
    Reads h5py file for given eigenvalue run. Plots the eigenfunctions.
    """
    fig, ax1 = plt.subplots()
    with h5py.File(file, 'r') as f:
        grid.S = f['Grid']['S'][()]
        r = grid.S
        for key in [x for x in f['Variables']['EvaluatedModes'].keys() if f'var{varnr}' in x]:
            mode = f['Variables']['EvaluatedModes'][key][()]
            mval = key.split("=",1)[1]
            ax1.plot(r, mode, label=f'm={mval}')

        ax1.legend(loc='best')
        ax1.set_title(f'var {varnr + 1}, file={file}')

    
def readh5(file, varnr = 0):
    """
    Reads h5py file for given eigenvalue run. Plots the eigenfunctions.
    """
    modes = {}
    with h5py.File(file, 'r') as f:
        grid.S = f['Grid']['S'][()]
        r = grid.S
        for key in [x for x in f['Variables']['EvaluatedModes'].keys() if f'var{varnr}' in x]:
            mode = f['Variables']['EvaluatedModes'][key][()]
            mval = key.split("=",1)[1]
            modes[f'm={mval}'] = mode
    return modes
    