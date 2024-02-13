import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import os
from datetime import datetime
import time
import random
import f90nml as f90
from pathlib import Path
from shutil import copy2

from VenusMHDpy import SATIRE2SFL
from VenusMHDpy import Equilibriumh5
from VenusMHDpy import Stability
from VenusMHDpy import VMECInput
from VenusMHDpy.library import find_nearest

class testScan(object):
    def __init__(self, runid = None, inputfile = None, inputpath = None):
        profiles = ['rho_t', 'rho_p', 'pt', 'rot']
        self.params = {'profile':'rho_t', 'mach':0.0, 'beta0':0.05, 'rmaj':10., 'a':1., 'b0':1., 'n':-1, 'rationalm':1, 'sidebands':5,
                        'd':0., 'el':0., 'tr':0., 'mpol':15, 'ntor':0, 'rstep':0.5, 'drstep':0.15, 
                        'n0':1, 'nu_n':2, 'qr':1., 'rs':0.3, 'q0':0.938, 'qs':None, 'nu_q':2.}
        
        self.initialisers = {'run_vmec': True, 'run_venus': True, 'toplot': True, 'peakedness': True, 
                        'ev_guess_type': 'last_ev', 'mode_type':'KH'}
        # EV_guess_type: determines what type of eigenvalue guessing system is used. 
        #'last_ev': last EV used as guess for next one, 'polynom_ev': polynomial fit used as guess for next one

        ### Define path to input file. Default is lunamhd/Input/default.in
        self.inputfile = inputfile
        if inputfile is None:
            self.inputfile = 'default.in'
        self.inputpath = inputpath
        if self.inputpath is None:
            #self.inputpath = f'/users/cs2427/lunamhd/Input/{self.inputfile}'
            self.inputpath = f'/home/csch/VENUS-linux/lunamhd/Input/{self.inputfile}'
        else:
            self.inputpath = self.inputpath + '/' + self.inputfile
        ### Define outpath
        if self['mode_type'] == 'KH':
            #self.outpath = Path('/users/cs2427/scratch/lunamhd-data/KH') # for running on viking
            self.outpath = Path('/home/csch/VENUS-linux/lunamhd/Output/KH/')
        elif self['mode_type'] == 'IK':
            #self.outpath = Path('/users/cs2427/scratch/lunamhd-data/IK') # for running on viking
            self.outpath = Path('/home/csch/VENUS-linux/lunamhd/Output/IK/')
        
        ### Essentially the contents of init_scan
        self.scans = {}
        self._readinput()
        if runid is None:
            runid = self._build_run_name(self.scanorder)
        self.runid = runid

    def __getitem__(self, key):
        if key in self.params:
            return self.params[key]
        elif key in self.initialisers:
            return self.initialisers[key]
        else:
            print(f"ERROR: {key} not found")
    
    ###### READ INPUT AND INITIALISE ######
    def _make_param_vals(self, nodes, res, descending = True):
        # Note: nodes should be provided in increasing order to avoid having negative paramvals (unless this is intended)
        # No error message for above since n values are negative
        steps = [nodes[0]]
        step_nrs = [1]
        for i in range(len(nodes)-1):
            step = (nodes[i] - nodes[i+1])/res[i]
            steps.append(step)
            step_nrs.append(res[i])
        deltas = np.repeat(steps, step_nrs)
        paramvals = np.cumsum(deltas)
        if descending:
            paramvals = np.flip(paramvals)
        return paramvals

    def _read_scanparams(self):
        # Scanorder: list of scan parameters with order retained (it matters which is 1st)
        # Scanparams: dict of scan parameters
        self.scanparams = {}
        self.scanorder = []

        scanparam_info = f90.read(self.inputpath)
        for key in [x for x in scanparam_info.keys() if 'scanparam_' in x]:
            scanparam_info = scanparam_info[key]
            paramkey = scanparam_info['varname']
            self.scanorder.append(paramkey)
            paramvals = self._make_param_vals(nodes = scanparam_info['nodes'], res = scanparam_info['resolution'], descending = scanparam_info['descending'])
            self.scanparams[paramkey] = paramvals

    def _readinput(self):
        # Reads input file.
        self._read_scanparams()

        inputs = f90.read(self.inputpath)
        inputs = inputs['config_nml']

        for key in inputs:
            if key not in ['scanparams','runid']:
                if key in self.params:
                    self.params[key] = inputs[key]
                elif key in self.initialisers:
                    self.initialisers[key] = inputs[key]
                else:
                    print(f"ERROR: {key} is invalid input, please review.")