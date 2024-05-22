# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:10:38 2023

@author: celin
"""
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
from copy import deepcopy

from VenusMHDpy import SATIRE2SFL
from VenusMHDpy import Equilibriumh5
from VenusMHDpy import Stability
from VenusMHDpy import VMECInput
from VenusMHDpy.library import find_nearest

class lunaScan(object):
    def __init__(self, runid = None, inputfile = None, inputpath = None):
        profiles = ['rho_t', 'rho_p', 'pt', 'rot']
        self.params = {'profile':'rho_t', 'mach':0.0, 'omegahat':0.0, 'beta0':0.05, 'rmaj':10., 'a':1., 'b0':1., 'n':-1, 'rationalm':1, 'sidebands':5, 'init_evguess':1E-1,
                        'd':0., 'el':0., 'tr':0., 'mpol':15, 'ntor':0, 'rstep':0.5, 'drstep':0.15, 
                        'n0':1, 'nu_n':2, 'qr':1., 'rs':0.3, 'q0':0.938, 'qs':None, 'nu_q':2.,
                            'grid_res':100}
        
        self.initialisers = {'run_vmec': True, 'run_venus': True, 'toplot': True, 'peakedness': True, 
                        'ev_guess_type': 'last_ev', 'mode_type':'KH', 'vmec_ver':8.5, 'splines':['q','p']}
        # EV_guess_type: determines what type of eigenvalue guessing system is used. 
        #'last_ev': last EV used as guess for next one, 'polynom_ev': polynomial fit used as guess for next one

        ### Default path stuff
        # self.inputpath_root = Path('/users/cs2427/lunamhd')
        self.inputpath_root = Path('/home/csch/VENUS-linux/lunamhd')
        # self.outputpath_root = Path('/users/cs2427/scratch/lunamhd-data')
        self.outputpath_root = Path('/home/csch/VENUS-linux/lunamhd/Output')

        ### Define path to input file. Default is lunamhd/Input/default.in
        self.inputfile = inputfile
        if inputfile is None:
            self.inputfile = 'default.in'
        self.inputpath = inputpath
        if self.inputpath is None:
            self.inputpath = Path(self.inputpath_root / f'Input/{self.inputfile}')
        else:
            self.inputpath = Path(self.inputpath + '/' + self.inputfile)
        
        ### Essentially the contents of init_scan
        self.scans = {}
        self._readinput()
        if runid is None:
            runid = self._build_run_name(self.scanorder)
        self.runid = runid

        ### Define outpath
        if self['mode_type'] == 'KH':
            self.outpath = Path(self.outputpath_root / 'KH')
        elif self['mode_type'] == 'IK':
            self.outpath = Path(self.outputpath_root / 'IK')

    def __getitem__(self, key):
        if key in self.params:
            return self.params[key]
        elif key in self.initialisers:
            return self.initialisers[key]
        else:
            print(f"ERROR: {key} not found")
    
    ###### READ INPUT AND INITIALISE ######
    def _make_param_vals(self, nodes, nsteps, descending = True):
        # Note: nodes should be provided in increasing order to avoid having negative paramvals (unless this is intended)
        # No error message for above since n values are negative
        steps = [nodes[0]]
        step_nrs = [1]
        if len(nodes) > 2:
            for i in range(len(nodes)-1):
                step = (nodes[i+1] - nodes[i])/nsteps[i]
                steps.append(step)
                step_nrs.append(nsteps[i])
        else:
            step = (nodes[-1] - nodes[0])/nsteps
            steps.append(step)
            step_nrs.append(nsteps)
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
            key_info = scanparam_info[key]
            paramkey = key_info['paramkey']
            self.scanorder.append(paramkey)
            if key_info['vals'] is None:
                paramvals = self._make_param_vals(nodes = key_info['nodes'], nsteps = key_info['nsteps'], descending = key_info['descending'])
                self.scanparams[paramkey] = paramvals
            else:
                paramvals = key_info['vals']
                self.scanparams[paramkey] = paramvals
        self.scandim = len(self.scanorder)

    def _readinput(self):
        # Reads input file.
        self._read_scanparams()

        inputs = f90.read(self.inputpath)
        inits = inputs['initialisers']
        fixed_params = inputs['fixed_params']

        for key in fixed_params:
            if key not in ['runid']:
                if key in self.params:
                    self.params[key] = fixed_params[key]
                else:
                    print(f"ERROR: {key} is invalid input, please review.")

        for key in inits:
            if key in self.initialisers:
                self.initialisers[key] = inits[key]
            else:
                print(f"ERROR: {key} is invalid initialiser, please review.")
        
        if self.initialisers['splines'] is None:
            self.initialisers['splines'] = []

    def _make_scan_list(self):
        # Makes the list of e.g. [{'beta':1,'delq':0.1},{'beta':1,'delq':0.2}]
        scandim = self.scandim
        
        def loop(n = scandim - 1, scandim = scandim, scanvars = {}, scans = []):
            if n == 0:
                return [{}]
            else:
                scanparam = self.scanorder[scandim - n]
                for val in self.scanparams[scanparam]:
                    scanvars[scanparam] = val
                    if n > 1:
                        loop(n = n - 1, scandim = scandim, scanvars = scanvars)
                    else:
                        scans.append(scanvars.copy())
            return scans
        return loop()

    def _build_scan_dir(self, scan):
        # Makes subdirectory for a given scan, e.g. Output/run/beta=1/delq=0.1/
        # Only relevant for integrated ND scans, not parallel ND scans (which function
        # as several 1D scans)
        run_saveloc = Path(f"{self.outpath}/{self.runid}")
        scandir = Path.cwd() / run_saveloc
        for param in scan:
            scandir = scandir / f'{param}_{scan[param]:.4f}'     
        return scandir

    def _build_scan_id(self, scan):
        # Makes scan ID in format beta_1.000-delq_0.100
        if scan:
            scanid = ""
            for paramkey, paramval in scan.items():
                        scanid += f"{paramkey}_{paramval:.4f}-"
        else:
            scanid = self.runid
        return scanid
    
    def _make_scan_inputs(self):
        # Makes scan input files and writes them to the appropriate subdirs
        # Also makes a self.list of the scan IDs + '\n'
        # Also makes a self.list of the subdir paths + '\n'
        self.scan_subdirs = []
        self.scan_ids = []
        for scan in self.scans:
            scan_subdir = self._build_scan_dir(scan)
            self.scan_subdirs.append(str(scan_subdir) + '\n')
            self.scan_ids.append(str(self._build_scan_id(scan) + '\n'))
            if not(Path.exists(scan_subdir)):
                scan_subdir.mkdir(parents=True)
            else:
                print("WARNING: Saving output to an existing file directory. Label is not unique?")
            old_inputs = f90.read(self.inputpath)
            for paramkey, paramval in scan.items(): # modify scanparams into their single values
                old_inputs['fixed_params'][paramkey] = paramval
            new_inputs = f90.Namelist()
            new_inputs['fixed_params'] = old_inputs['fixed_params']
            new_inputs['initialisers'] = old_inputs['initialisers']
            new_inputs['scanparam_0'] = old_inputs['scanparam_0']
            with open(scan_subdir / self.inputfile, 'w') as f:
                f90.write(new_inputs, f, force=True)
            
    def _write_scan_dirs(self):
        # Writes the scan subdir paths to a file 'scan_subdirs.txt'
        # Writes the scan IDs to a file 'scan_ids.txt'
        run_saveloc = Path(f"{self.outpath}/{self.runid}")
        inputs_list_file = run_saveloc / f'scan_subdirs.txt'
        with open(inputs_list_file, 'w') as f:
            f.writelines(self.scan_subdirs)
        scanids_list_file = run_saveloc / f'scan_ids.txt'
        with open(scanids_list_file, 'w') as f:
            f.writelines(self.scan_ids)

    def init_run(self):
        # Initialises a run.
        self.scans = self._make_scan_list()
        self._make_scan_inputs()
        self._write_scan_dirs()
        return
            
    ###### VMEC ######
    def _peakedness1(self, r, y):
        # Simplified peakedness calculation
        
        dlprof = len(y)-99 # to force the length of y to match r
        dldprof = len(y)-100
        
        prof = y[dlprof:] # profile is equal to array of values given, skip y(a) since prof needs to be same length as dprof
        dprof = [i/j for i, j in zip(np.diff(y[dldprof:]), np.diff(r))]
        
        p1 = spi.trapezoid(r[1:]*dprof/prof)
        
        return p1
    
    def _peakedness2(self, r, y, xi):
        # Takes an array as y input. 
        # Potential improvement: profiles are approximated as polynomials --> use the actual polynomial derivative and stuff? But sometimes they're spline
        
        prof = y[1:] # profile is equal to array of values given, skip y(a) since prof needs to be same length as dprof
        dprof = [i/j for i, j in zip(np.diff(y), np.diff(r))]
        
        print(len(xi))
        print(len(r[1:]))
        print(len(dprof))
        p2 = spi.trapezoid(xi[1:]*r[1:]*dprof/prof)/max(xi)    
        
        return p2

    def _build_VMEC_profile(self, param):
        # able to pick between spline or polynomial
        # special conditions for e.g. density(?) needing to have nonzero bit
        # can tell it which parameter and it produces the right variable?
        return

    def _buildVMEC(self, idx = 0):
              
        #Read the default input file
        C = VMECInput.ReadInputVMEC(self.inputpath_root / 'VMEC/input/input.Default') 
        
        #Modify some grid and control parameters, these get written to VMEC input
        #======================================================================
        C.Grid.MPOL = self.params['mpol']    #Number of poloidal modes used
        C.Grid.NTOR = self.params['ntor']    #Set 2D	
        C.Grid.LASYM = 'F'	#Weather or not to violate stellarator symmetry.
        C.Grid.NZETA = 16   #Number ot toroidal planes. Important in free boundary calculations
        C.Grid.LRFP = 'F'   #Weather to use toroidal (F) or poloidal (T) normalized flux as radial variable. Note that if poloidal flux is used, 'q' needs to be provided instead of 'iota'.
        C.FreeB.LFREEB = 'F'    #Set the simulation to be fixed boundary
        #======================================================================
        
        #Boundary
        #======================================================================
        mu0 = 4.*np.pi*1.0E-07
        e = 1.60217663E-19
        r0 = self['a']
        R0 = self['rmaj']
        B0  = self['b0']
        D = self['d']  #Shafranov Shift
        El = self['el'] #Elongation
        Tr = self['tr'] #Triangularity limit 0.08
        
        n = [0,0,0] #n and m are linked to the mode nb for R and Z at the boundaries 
        m = [0,1,2]
        RBC = [R0,r0-El,Tr]
        ZBS = [0.,r0+El,-Tr]
        RAXIS = [RBC[0]]
        ZAXIS = [0.]
        
        C.Boundary.PHIEDGE = np.pi*r0**2.*B0
        C.Boundary.RAXIS = RAXIS
        C.Boundary.ZAXIS = ZAXIS
        C.Boundary.RBCn  = n
        C.Boundary.RBCm  = m
        C.Boundary.RBC   = RBC
        C.Boundary.ZBSn  = n
        C.Boundary.ZBSm  = m
        C.Boundary.ZBS   = ZBS
        #======================================================================
        
        
        #Grid
        #==========================
        s2 = np.linspace(0.,1.,100)
        s  = np.sqrt(s2)
        #==========================
        
        
        ### PARAMETER PROFILES ###
        #######################################################################
        
        ### Density: not a VMEC input but useful so we can set pressure wrt density.
        ### PT & rot: output AT and AH will be same as input but for the 
        # computation VMEC will normalize the profiles as T --> T/T0, 
        # Omega --> Omega/Omega0
        ### Pressure: care with P vs PVMEC

        mach = self['mach']
        beta0 = self['beta0']
        n0 = self['n0']
        nu_n = self['nu_n']

        if self['mode_type'] == 'KH':
            rstep = self['rstep']
            drstep = self['drstep']

            if self['profile'] != 'rot':
                Omega = (1.-s**6)
                
                if self['profile'] in ['rho_t', 'rho_p']:
                    n_ = .5*n0*(1 + np.tanh((rstep**2 - s2)/drstep**2))
                    
                    if self['profile'] == 'rho_t':
                        n_ = n_ + 0.05
                        T = .5*(1 - np.tanh((rstep**2 - s2)/drstep**2)) + 0.05
                        P = beta0*B0**2/(2*mu0*n0*T[0]) # should be constant
                        #P = beta0*B0**2/(2*mu0*n0*T[0])*(1-s2**3)
                        
                    elif self['profile'] == 'rho_p':
                        T = np.ones_like(s)
                        P = beta0*B0**2*n_*T/(2*mu0*n0) # stepped like n_ 

                elif self['profile'] == 'PT':
                    n_ = n0*np.ones_like(s)
                    T = .5*(1 + np.tanh((rstep**2 - s2)/drstep**2)) + 0.05
                    P0 = beta0*B0**2*n_*T/(2*mu0*n0)
                    P = .5*P0*(1 + np.tanh((rstep**2 - s2)/drstep**2))
                
            elif self['profile'] == 'rot':
                Omega = .5*(1 + np.tanh((rstep**2 - s2)/drstep**2))
                n_ = n0*(1.-s**nu_n)
                T = n0*(1.-s**6) + 0.01
                P = beta0*B0**2*n_*T/(2*mu0*n0)  
        
        elif self['mode_type'] == 'IK':
            n_ = n0*1.-s**nu_n
            # n_ = np.ones_like(s)
            T = np.ones_like(s) # makes density equal to the pressure
            # T = 1-s2+0.05
            P = beta0*B0**2.*(1-s2)/(2.*mu0)
            Omega = np.ones_like(s)
        
        ### ROTATION
        Omega = Omega/Omega[0]
        self.Omega = Omega

        if 'o' in self['splines']:
            print('=== SPLINING FLOW ===')
            C.Flow.PH_TYPE = "'cubic_spline'"
            AH_AUX_S = s2
            AH_AUX_F = Omega
            AH = [0]
        else:
            AH_AUX_S = [0]
            AH_AUX_F = [0]
            AH = np.polyfit(s2,Omega,11)[::-1]

        C.Flow.AH_AUX_S = AH_AUX_S
        C.Flow.AH_AUX_F = AH_AUX_F
        C.Flow.AH = AH # SET FLOW PROFILE
        C.Flow.bcrit = mach # SET FLOW MAGNITUDE
        
        ### PRESSURE (order matters to get normalized Omega)
        PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)

        if 'p' in self['splines']:
            print('=== SPLINING PRESSURE ===')
            C.Pressure.PMASS_TYPE = "'cubic_spline'"
            C.Pressure.AM_AUX_S = s2
            C.Pressure.AM_AUX_F = PVMEC
        else:
            AM = np.polyfit(s2,PVMEC,11)[::-1]
            C.Pressure.AM = AM 
            C.Pressure.PRES_SCALE = 1.
        
        # ### TEMPERATURE
        T = T/T[0]

        if 't' in self['splines']:
            print('=== SPLINING TEMPERATURE ===')
            C.Flow.PT_TYPE = "'cubic_spline'"
            AT_AUX_S = s2
            AT_AUX_F = T
            AT = [0]
        else:
            AT_AUX_S = [0]
            AT_AUX_F = [0]
            AT = np.polyfit(s2,T,11)[::-1]

        C.Flow.AT_AUX_S = AT_AUX_S
        C.Flow.AT_AUX_F = AT_AUX_F
        C.Flow.AT = AT

        ### Q PROFILE
        #======================================================================
        C.Current.NCURR  = 0     #0 for rotal transform, 1 for toroidal current density
        #======================================================================
        qr = self['qr']
        rs = self['rs'] # set to 0 to get qs = 1, this is r where q = qr, equivalent to r1 in Tom's IK work
        q0 = self['q0']
        nu_q = self['nu_q']
        
        if rs == 0:
            qs = self['qs']
            q = q0+qs*s**nu_q
        else:
            qs = (qr-q0)/rs**(nu_q)
            q = 1 - (1-q0)*(1 - (s/rs)**nu_q)
            
        if 'q' in self['splines']:
            print('=== SPLINING Q-PROFILE ===')
            C.Current.PIOTA_TYPE = "'cubic_spline'"
            C.Current.AI_AUX_S = s2
            if C.Grid.LRFP == 'F':
                C.Current.AI_AUX_F = -1./q
            elif C.Grid.LRFP == 'T':
                C.Current.AI_AUX_F = -q
            else:
                print ('Insert a valid value for LRFP')
                exit()
        else:
            if C.Grid.LRFP == 'F':
                AI = np.polyfit(s2,-1./q,11)[::-1]
            elif C.Grid.LRFP == 'T':
                AI = np.polyfit(s2,-q,11)[::-1]
            else:
                print ('Insert a valid value for LRFP')
                exit()
            C.Current.AI = AI

        self.dico_vmec = {'ah':AH,'ah_aux_s':AH_AUX_S,'ah_aux_f':AH_AUX_F,'at':AT,'at_aux_s':AT_AUX_S,'at_aux_f':AT_AUX_F}
        
        #Change some control parameters
        #======================================================================
        """
        C.Control.PREC2D_THRESHOLD = 1.0E-13 
        C.Control.NITER_ARRAY = [1999, 3999, 3999, 3999, 8999, 8999, 8999, 8999, 8999, 25999, 39999, 99999, 129999]
        C.Control.NS_ARRAY    = [25, 73, 211, 321, 435, 449, 463, 475, 481, 483, 485, 487, 489]
        C.Control.FTOL_ARRAY  = [1.0e-09, 1.0e-09, 5.0e-10, 5.0e-10, 5.0e-10, 1.0e-10, 5.0e-11, 5.0e-11, 5.0e-11, 5.0e-11, 1.0e-11, 1.0e-11, 5.0E-12]
        """
        #======================================================================
        
        #Run VMEC Fixed boundary VMEC
        #======================================================================
        DIR_VMEC = f'{self.inputpath_root}/VMEC/' # BUGFIX: MAKE SURE THIS IS SET CORRECTLY
        Fout = 'input.'+f'{self.runid}_{idx}'
        if self['run_vmec']:
        	#Write the input file
            C.WriteInput(Fout)
        
        	#Run VMEC
            if self['vmec_ver'] == 8.5:
                vmec_ver = 'xvmec2000_flow_netcdf'
            elif self['vmec_ver'] == 9:
                vmec_ver = 'vmec_flow'
            else:
                print(f'ERROR: invalid VMEC version specified.')
            os.system(DIR_VMEC+f'./{vmec_ver} '+Fout)
            print(f"USING VMEC VERSION {vmec_ver}")

            #Create the folders if not existent
            input_f = Path(self.outpath / f'{self.runid}/VMEC/input')
            mercier_f = Path(self.outpath / f'{self.runid}/VMEC/mercier')
            jxbout_f = Path(self.outpath / f'{self.runid}/VMEC/jxbout')
            threed1_f = Path(self.outpath / f'{self.runid}/VMEC/threed1')
            wout_f = Path(self.outpath / f'{self.runid}/VMEC/wout')

            input_f.mkdir(parents=True, exist_ok=True)
            mercier_f.mkdir(parents=True, exist_ok=True)
            jxbout_f.mkdir(parents=True, exist_ok=True)
            threed1_f.mkdir(parents=True, exist_ok=True)
            wout_f.mkdir(parents=True, exist_ok=True)

            #Move the output files to their folders

            Path(f'input.{self.runid}_{idx}').rename(input_f / f'input.{self.runid}_{idx}')
            Path(f'wout_{self.runid}_{idx}.nc').rename(wout_f / f'wout_{self.runid}_{idx}.nc')
            Path(f'mercier.{self.runid}_{idx}').rename(mercier_f / f'mercier.{self.runid}_{idx}')
            Path(f'jxbout_{self.runid}_{idx}.nc').rename(jxbout_f / f'input.{self.runid}_{idx}')
            Path(f'threed1.{self.runid}_{idx}').rename(threed1_f / f'input.{self.runid}_{idx}')
            os.system('rm dcon_'+f'{self.runid}_{idx}'+'.txt')
        #======================================================================
        
    ###### VENUS ######
    def _runVENUS(self, EVguess = None, idx = 0):
        
        """
        Returns Gamma/OmegaA.
        
        Parameters:
            EV_guess - Initial guess for eigenvalue calculation
            idx - Index for scans. If idx = 0, default EV_guess is used. If
            idx > 0, previous eigenvalue calculated in scan is used.
            labelnr - Does not affect eigenvalue guess. Is used to produce plots
            for a specific single sweep within a scan over several parameter
            values by picking 1 specific VMEC input file.
        """
        
        #Read equilibrium from VMEC output file and transform it into SFL
        eq = SATIRE2SFL.SATIRE2SFL(woutfile = self.outpath / f'{self.runid}/VMEC/wout/wout_{self.runid}_{idx}.nc', dico_vmec = self.dico_vmec)
        #eq = SATIRE2SFL.SATIRE2SFL(woutfile = self.outpath / f'{self.runid}/VMEC/wout/wout_{self.runid}_{idx}.nc')
        eq.Writeh5(eqFile=f'eq.{self.runid}_{idx}.h5')
        os.system('mv '+f'eq.{self.runid}_{idx}.h5'+' eqFiles')
    	
    	#Create the stability object
        stab = Stability.Stability('IdealMHDFlow-Euler')
        eq.kappa = 0.
        
    	#Modify the default grid
    	#----------------------------------------------------------------------
        n = self['n'] # toroidal mode number, <0 because of how vars are expanded in n, m
        n = -1
        RationalM = self['rationalm']
        RationalM = 1
        Sidebands = self['sidebands']
        Sidebands = 5
        stab.grid.Mmin = RationalM-Sidebands
        stab.grid.Mmax = RationalM+Sidebands
        stab.grid.Ntheta = eq.R.shape[0]
    
        stab.grid.N = 100
        stab.grid.N = eq.R.shape[1] - 2
        stab.grid.bunching = False
        stab.grid.bunchingQValues = [1.0,1.1, 1.2]
        stab.grid.bunchingAmplitudes = [5.,5.,5.]
        stab.grid.bunchingSigma = [0.02,0.02,0.02]
        
    	#Build grid. If bunching with q values, then equilibrium quantities (radial grid s and safety factor) are required.
        stab.grid.BuildGrid(eq.s,eq.q)
    	#----------------------------------------------------------------------
        
    	#Normalize and build the equilibrium quantities in the new grid.
    	#----------------------------------------------------------------------
        eq.ChangeGrid(stab.grid.S)
        eq.Normalise()
        eq.BuildInGrid(stab.grid)

        if self['mode_type'] == 'IK':
            eq.Omega = -eq.Omega # KH doesn't run well if this is not set

        # Calculate Shafranov stuff
        LHS = eq.dFds*eq.g22/(eq.q*eq.R2)+eq.F*eq.dg22ds/(eq.q*eq.R2)-eq.F*eq.g22*eq.dqds/(eq.q**2*eq.R2)-eq.F*eq.g22*eq.dR2ds/(eq.q*eq.R2**2.)-eq.F/eq.q*(eq.dg12du/eq.R2-eq.g12*eq.dR2du/eq.R2**2.)
        RHS = -eq.q*eq.dFds-(eq.q*eq.R2/eq.F)*(eq.dPds+(eq.R**2.-1.)*eq.P*eq.dUds)*np.exp(eq.U*(eq.R**2.-1.))
        shaf_diff0 = (RHS[0]-LHS[0])/max(abs(LHS[0]))
    	
        V0_Va = np.sqrt(eq.M02*eq.mu0*eq.P0)/eq.B0
        # Calculate Omegahat, rotation frequency as normalized in 2013 PPCF
        eps_a = self['a']/self['rmaj']
        betahat = 2*eq.mu0*eq.P0/eq.B0**2/eps_a**2
        Omegahat = np.sqrt(eq.M02*2*eq.mu0*eq.P0)/(eq.B0*eps_a)
        self.params['omegahat'] = Omegahat
        # Print equilibrium quantities
        print ('Parameters at the magnetic axis:')
        print ('   M0    = %.5f'%(np.sqrt(eq.M02)))
        print ('   v0/vA = %.5f'%(V0_Va))
        print ('   B0    = %.5f / %.5f [T]'%(eq.B0, self['b0']))
        print ('   R0    = %.5f / %.5f [m]'%(eq.R0, self['rmaj']))
        print ('   P0    = %.5f / %.5f [Pa]'%(eq.P0, (self['beta0']*self['b0']**2/(2.*eq.mu0)))) #P = beta0*B0**2.*n_*T/(2.*mu0*n0)
        print ('   beta0 = %.5f / %.5f %%'%(2.*eq.mu0*eq.P0/eq.B0**2., self['beta0']))
    	#----------------------------------------------------------------------
        
        if True:
    		# Discretize the Operators.
    		#------------------------------------------------------------------
            stab.Discretize(eq, n)
    		#------------------------------------------------------------------
    
    		# Solve
    		#------------------------------------------------------------------
            t0 = time.time()
            if EVguess == None:
                idx_rstep = find_nearest(stab.grid.S, self['rstep'])
                #EV_guess = 1.0E-1 + (1.0j)*abs(n)*eq.Omega[idx_rstep] # want to re-implement this
                EVguess = self['init_evguess']
                EVguess = EVguess + (1.0j)*abs(n)*eq.Omega[0]
            #elif EV_guess == 'bad': #EV_guess.real < 1.0E-07 # an attempt at correcting when the EV guesses get bad
                #idx_rstep = find_nearest(stab.grid.S, self.profParams['rstep'])
                #EV_guess = 1.0E-3 + (1.0j)*abs(n)*eq.Omega[idx_rstep]
                
            print ('EV guess: %.5E + i(%.5E)'%(EVguess.real,EVguess.imag))
            stab.Solve(EVguess,N_EV=1)
            stab.Saveh5(FileName=f'{self.scan_saveloc}/{self.runid}_{idx}', eq = eq) # save eigenfunctions + eigenvalues + run info
            print ('Solution time: %.4E with N = %i' %(time.time()-t0, stab.grid.N))
    		#------------------------------------------------------------------
    
            EV = max(stab.vals)
            
            # Calculate peakedness
            # s2 = np.linspace(0.,1.,100) # gridblz
            # if self.initialisers['Peakedness']:
            #     if self.params['profile'] in ['rho_t', 'rho_p']:
            #         #P = beta0*B0**2.*n_*T/(2.*mu0*n0)
            #         mu0 = 4.*np.pi*1.0E-07
            #         dens = eq.P*2*mu0/(eq.beta0*eq.B0**2*eq.T) # normalised density
            #         pkedness = self._peakedness1(s2, dens)                        
            #     elif self.params['profile']=='PT':
            #         #pkedness = self._peakedness2(s2, self.Tpoly, xi)
            #         pkedness = self._peakedness1(s2, eq.T)
            #     elif self.params['profile']=='rot':
            #         #pkedness = self._peakedness2(s2, self.rotpoly, xi)
            #         pkedness = self._peakedness1(s2, eq.Omega)
            #     else:
            #         print('PROFILE IS NOT SET TO ONE OF THE ESTABLISHED STEPPED PROFILES.')
            
    		
            print ('Most unstable eigenvalue')
            print ('(Gamma/OmegaA) = %.5E + i(%.5E)'%(EV.real,EV.imag))

        output = {'EV':EV, 'v0_va':V0_Va, 'EVguess':EVguess, 'EF_file':f'{self.scan_saveloc}/{self.runid}_{idx}.h5', 'params':self.params.copy(), 'profile':self['profile'], 'shaf_diff0':shaf_diff0}
        
        return output.copy()

    def firstscan(self, scanid = 'default', scanparam = None):
        if scanparam is None:
            try:
                scanparam = self.scanorder[0]
            except:
                print("ERROR: scanparams is likely empty, check input file has read correctly")
        
        output1d = {}
        for vidx, val in enumerate(self.scanparams[scanparam]):
            self.params[scanparam] = val # key needs to be the same in params and scanparams
            # Run VMEC
            if self['run_vmec']:
                self._buildVMEC(idx = vidx) # sets runid inside of this function
                eq = SATIRE2SFL.SATIRE2SFL(woutfile = self.outpath / f'{self.runid}/VMEC/wout/wout_{self.runid}_{vidx}.nc', dico_vmec = self.dico_vmec)
                #eq = SATIRE2SFL.SATIRE2SFL(woutfile = self.outpath / f'{self.runid}/VMEC/wout/wout_{self.runid}_{vidx}.nc')

            if self['run_venus']:
                # Set EV guess and calculate the growth rate
                if self['ev_guess_type'] == 'last_ev':
                    if vidx <= 1:
                        output1d[f'{scanid}_{vidx}'] = self._runVENUS(EVguess = None, idx = vidx)
                    else:
                        lastrun = output1d[f'{scanid}_{vidx-1}']
                        EVguess = lastrun['EV']
                        output1d[f'{scanid}_{vidx}'] = self._runVENUS(EVguess = EVguess, idx = vidx)
                elif self['ev_guess_type'] == 'polynom_ev':
                    ws = [] # BUGFIX: this is probably gonna be a problem for re-running scans halfway through, going to need to find a way to read results as they're being made
                    if vidx <= 2:
                        output1d[f'{scanid}_{vidx}'] = self._runVENUS(EVguess = None, idx = vidx)
                    else:
                        polycoeff = vidx - 2
                        if polycoeff > 10: # changing these numbers can help improve fits sometimes
                            polycoeff = 10
                        ws = []
                        for i in range(vidx):
                            ws.append(output1d[f'{scanid}_{i}']['EV'])
                        guessReal = np.polyfit(np.asarray(self.scanparams[scanparam][:vidx]),np.asarray([i.real for i in ws]),polycoeff)
                        guessImag = np.polyfit(np.asarray(self.scanparams[scanparam][:vidx]),np.asarray([i.imag for i in ws]),polycoeff)
                        
                        EVguess = np.polyval(guessReal,val)*3 + 1j*np.polyval(guessImag,val)
                        EVguess += 1j*EVguess.imag*1E-3 # want to be slightly larger than the correct EV
                        output1d[f'{scanid}_{vidx}'] = self._runVENUS(EVguess = EVguess, idx = vidx)
                        ws.append(output1d[f'{scanid}_{vidx}']['EV'])
        
        return output1d.copy()

    def run(self, scan_saveloc = None):
        self.scan_saveloc = scan_saveloc
        if self.scans: # integrated ND scan
            for scan in self.scans:
                scanoutput = self.firstscan(scanid = self.runid)
                self._save_scan(scanid = self.runid, scanoutput = scanoutput, scan_saveloc = self.scan_saveloc)
        else: # 1D (or parallel ND) scan
            scanoutput = self.firstscan(scanid = self.runid)
            self._save_scan(scanid = self.runid, scanoutput = scanoutput, scan_saveloc = self.scan_saveloc)        
        return
    
    ###### BUILD OUTPUT FILE ######
    def _get_scan_info(self, scanid = None):
        scan_info = {}
        scan_info['scanparams'] = self.scanparams        
        scan_info['timestamp'] = datetime.now().strftime("%d-%m-%y_%H;%M")
        return scan_info.copy()

    def _build_run_name(self, scankeys):
            # Generates a run name e.g. rot_Omega_beta_delq
            profile = self['profile']
            out_filename = f"{profile}"
            for key in scankeys:
                out_filename += f"_{key}"     
            return out_filename

    def _save_scan(self, scanoutput, scanid = None, scan = None, scan_saveloc = None):
        # Saves a single scan to its corresponding subdirectory
        scan_info = self._get_scan_info(scanid = self.runid)
        if scan is None:
            if self.scan_saveloc is None:
                fOut = self.outpath / f'{self.runid}.npz' # self.scanid.npz?
            else:
                fOut = Path(self.scan_saveloc) / f'{self.runid}.npz' # self.scanid.npz?
        else:
            fOut = self._build_scan_dir(scan) / f'{self.runid}.npz' # self.scanid.npz?
        if not(Path.exists(self.outpath)):
            self.outpath.mkdir(parents=True)
            np.savez(fOut, data = scanoutput.copy(), info = scan_info.copy())
        else:
            np.savez(fOut, data = scanoutput.copy(), info = scan_info.copy()) 
        return

    def save_run(self, runid = None):
        if runid is None:
            runid = self.runid
        run_saveloc = Path(f"{self.outpath}/{self.runid}")
            
        rundata = {}
        runinfo = {}
        scans = deepcopy(self.scans)
        if self.scans:
            print(f'self.scans = {self.scans}')
            for n, scan in enumerate(self.scans):
                scandir = self._build_scan_dir(scan)
                scanid = f'{runid}_{n}'
                try:
                    scanfile = sorted(scandir.glob('*.npz'))[0] # picks out the first .npz file
                    raw_scan = np.load(scanfile, allow_pickle = True)
                    rundata[scanid] = raw_scan['data'].item()
                except:
                    print(f'Scan {scan} output file not found. May not have converged.')
                    scans.remove(scan) # remove the skipped scan from the scanlist for the final output file
                    continue
                print(scan)
        else:
            runfile = sorted(run_saveloc.glob('*.npz'))[0]
            raw_scan = np.load(runfile, allow_pickle = True)
            rundata[runid] = raw_scan['data'].item()
            
        runinfo['scanparams'] = self.scanparams
        fixed_params = ['profile','drstep','mach','q0','beta0','rationalm','n', 'rmaj']
        runinfo['fixedparams'] = {}
        for key in fixed_params:
            if key not in self.scanparams:
                runinfo['fixedparams'][key] = self[key]
        runinfo['scanorder'] = self.scanorder
        runinfo['scans'] = self.scans
        runinfo['timestamp'] = datetime.now().strftime("%d-%m-%y_%H:%M")
        runinfo['runid'] = self.runid
        
        fOut = f"{run_saveloc}/{self.runid}.npz"
        
        np.savez(fOut, data = rundata, info = runinfo)
        # add any new info to the final inputfile
        inputfile = Path.cwd() / self.inputpath
        inputfile_nml = f90.read(inputfile)
        inputfile_nml['info']['runid'] = self.runid
        inputfile_nml['info']['datetime'] = datetime.now().strftime("%d-%m-%y_%H:%M")
        with open(run_saveloc / self.inputfile, 'w') as f:
            f90.write(inputfile_nml, f, force=True)
        #copy2(inputfile, run_saveloc)
        
        return # can save old runs if same input file (needs same scan parameters) and runid is provided using init_run

    def local_run(self):
        self.init_run()
        self.run(self.outpath / self.runid)
        self.save_run()
                    

        
    
    