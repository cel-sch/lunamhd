# -*- coding: utf-8 -*-
"""
Created on Tue May 30 18:59:04 2023

@author: celin
"""

import h5py
import numpy as np
from matplotlib.pyplot import subplots, show, axes, rcParams, figure
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

    
def readh5_EF(file, varnr = 0):
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
    
def readh5_profs(file, show = True):
    with h5py.File(file, 'r') as f:
        mu0 = 4.*np.pi*1.0E-07

        # Unpack variables
        B0 = f['normalisation']['B0'][()]
        P0 = f['normalisation']['B0'][()]
        R0 = f['normalisation']['B0'][()]
        M02 = f['normalisation']['M02'][()]

        s = f['profiles']['s'][()]
        T = f['profiles']['T'][()]
        q = f['profiles']['q'][()]
        P = f['profiles']['P'][()]
        Prot = f['profiles']['Prot'][()]
        Omega = f['profiles']['Omega'][()]
        U = f['profiles']['U'][()]
        rho = f['profiles']['rho'][()]

        R = f['geometry']['R'][()]
        Z = f['geometry']['Z'][()]

    def ploth5_profs(show=True):
        rcParams['axes.grid'] = True

        fig = figure()
        ax1 = fig.add_subplot(2,4,1)
        ax1.plot(s, 2*B0**2*T/(P0*mu0))
        ax1.set_ylabel(r'$T/T_0$')
        ax1.set_xlabel('s')

        ax2 = fig.add_subplot(2,4,2)
        ax2.plot(s, q)
        ax2.set_ylabel(r'q')
        ax2.set_xlabel('s')

        ax3 = fig.add_subplot(2,4,3)
        ax3.plot(s, B0**2*P/mu0, label=r'$P(s)$')
        ax3.plot(s, B0**2*Prot[0]/mu0, label=r'$P(s)e^{U(R^2-R_0^2)}$')
        ax3.legend()
        ax3.tick_params('x', labelbottom=False)
        ax3.set_ylabel(r'$\bar{P}$ [Pa]')
        ax3.set_xlabel('s')

        ax4 = fig.add_subplot(2,4,4)
        ax4.plot(R[:,-1]*R0, Z[:,-1]*R0)
        ax4.plot(R[0,0]*R0, Z[0,0]*R0,'+')
        ax4.set_xlabel(r'R [m]')
        ax4.set_ylabel(r'Z [m]')
        ax4.set_aspect('equal')

        ax5 = fig.add_subplot(2,4,5)
        ax5.plot(s, U/R0**2.)
        ax5.set_ylabel(r'$U$ $[m^{-2}]$')
        ax5.set_xlabel('s')

        ax6 = fig.add_subplot(2,4,6)
        ax6.plot(s, rho)
        ax6.set_ylabel(r'$\bar{\rho}/\rho_0$')
        ax6.set_xlabel('s')

        ax7 = fig.add_subplot(2,4,7)
        if M02 == 0:
            ax7.plot(s,np.zeros_like(s))
        else:
            ax7.plot(s,Omega*B0/(M02*mu0*P0)**0.5)
        ax7.set_ylabel(r'$Ω/Ω_0$')
        ax7.set_xlabel('s')

        if show:
            show()

    ploth5_profs(show=show)

    return

