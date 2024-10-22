# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:10:38 2023

@author: celin
"""
import h5py
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
        profiles = ['rot_0', 'all_step', 'all_step_tvo']
        self.params = {'profile':'rho_t', 'mach':0.0, 'omegahat':0.0, 'rmaj':10., 'a':1., 'b0':1., 'n':-1, 'rationalm':1, 'sidebands':5, 'init_evguess':1E-1,
                       'omega_avg':1, 'omega_step':0.9, 'omega0':0, 'omega1':0,
                       'rho_avg':1, 'rho_step':0, 'rho0':0, 'rho1':0,
                       'beta_avg':0.05, 'beta_step':0, 'beta0':0, 'beta1':0,
                        'd':0., 'el':0., 'tr':0., 'mpol':15, 'ntor':0, 'rstep':0.5, 'drstep':0.15, 
                        'n0':1, 'nu_n':2, 'qr':1., 'rs':0.3, 'delq':0.938, 'qs':None, 'nu_q':2.,
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
            self.inputfile = 'default_KH.in'
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
        mode_type = f"{self['mode_type']}"
        self.outpath = Path(self.outputpath_root / mode_type)

        ### Copy input file to output location
        self._copyinput()

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
            if paramkey not in self.params:
                print('ERROR: specified scan parameter does not correspond to any of the possible input parameters.')
                return
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

    def _copyinput(self):
        run_saveloc = Path(f"{self.outpath}/{self.runid}")
        run_saveloc.mkdir(parents=True, exist_ok=True)
        # add any new info to the final inputfile
        inputfile = Path.cwd() / self.inputpath
        inputfile_nml = f90.read(inputfile)
        inputfile_nml['info']['runid'] = self.runid
        inputfile_nml['info']['datetime'] = datetime.now().strftime("%d-%m-%y_%H:%M")
        with open(run_saveloc / self.inputfile, 'w') as f:
            f90.write(inputfile_nml, f, force=True)

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
            new_inputfile = scan_subdir / self.inputfile
            if new_inputfile.is_file():
                pass
            else:
                with open(scan_subdir / self.inputfile, 'w') as f: # is this overwriting the copied main input file? yes
                    f90.write(new_inputs, f)
            
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
    def _peakedness(self, s, y, dyds, xi, figidx):
        y = np.asarray(y)
        dyds = np.asarray(dyds)
        xi = np.asarray(xi)
        xi_ = np.asarray(xi)
        xi[1:] = np.asarray([i/j for (i, j) in zip(xi[1:], s[1:])]) # remove *s dependence that's probably in var1, skip s=0 point
        xi = abs(np.asarray([i/max(abs(xi)) for i in xi])) # successfully normalises xi to 1
        p2 = abs(spi.simps(xi*dyds, x=s))/np.mean(y) # with y_avg normalisation
        p = abs(spi.simps(xi*dyds, x=s)) # without y_avg normalisation

        print(f"integral dy/ds: {spi.simps(dyds, x=s)}")
        print(f"peakedness: {p}")

        if figidx == 0:
            plt.figure()
            plt.plot(s, xi)
            plt.savefig(f"{self.runid}_xi.png")
        return p, p2

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
        s2 = np.linspace(0.,1.,100) # can't change the grid size when splining?
        s  = np.sqrt(s2)
        #==========================
        
        
        ### PARAMETER PROFILES ###
        #######################################################################
        
        ### Density: not a VMEC input but useful so we can set pressure wrt density.
        ### PT & rot: output AT and AH will be same as input but for the 
        # computation VMEC will normalize the profiles as T --> T/T0, 
        # Omega --> Omega/Omega0
        ### Pressure: care with P vs PVMEC

        def step_profile(y_avg, y_step):
            # for a tanh profile of the form 0.5*B(1 + tanh[(r0^2-r^2)/dr0^2]) + A
            if y_step == 0:
                #prof = y_avg*(1.-s**6)
                prof = y_avg*np.ones_like(s)
            else:
                B = 2*y_step
                A = y_avg - y_step
                if self['profile'] == 'rot': # Omega1 = 0 case with all other profiles flat
                    A = 0
                    B = y_step

                prof = B*0.5*(1 + np.tanh((rstep**2 - s2)/drstep**2)) + A
                if any(i < 0 for i in prof):
                    print('ERROR: profile contains a negative value. Check y_avg and y_step are correctly specified (did you swap them?).')
            return prof
        
        def make_y0_y1(y_avg, y_step):
            y0 = y_avg + y_step
            y1 = y_avg - y_step
            return y0, y1

        mach = self['mach']
        beta0 = self['beta_avg']
        n0 = self['n0']
        nu_n = self['nu_n']

        if self['mode_type'] == 'KH':
            rstep = self['rstep']
            drstep = self['drstep']

            self.params['omega0'], self.params['omega1'] = make_y0_y1(self['omega_avg'], self['omega_step'])
            self.params['rho0'], self.params['rho1'] = make_y0_y1(self['rho_avg'], self['rho_step'])
            self.params['beta0'], self.params['beta1'] = make_y0_y1(self['beta_avg'], self['beta_step']) # if beta_step=0, beta0=beta_avg

            if self['profile'] in ['all_step', 'all_step_tvo']: # rhostep = 0 with stepped omega, pressure is fine but having flat density AND pressure is bad so farpfp 
                beta0 = self.params['beta0']

                Omega = step_profile(self['omega_avg'], self['omega_step'])
                n_ = step_profile(self['rho_avg'], self['rho_step']) + 0.01 # need +0.01 for 1/n_ in T

                if self['profile'] == 'all_step_tvo':
                    beta1 = beta0*(self['rho1']/self['rho0'])*(self['omega1']/self['omega0'])**4
                    self.params['beta_avg'] = (beta0 + beta1)/2
                    self.params['beta_step'] = (beta0 - beta1)/2
                    self.params['beta0'], self.params['beta1'] = make_y0_y1(self['beta_avg'], self['beta_step'])

                beta = step_profile(self['beta_avg'], self['beta_step'])    
                P = beta*B0**2/(2*mu0) 
                T = 2*mu0*n_[0]*P/(beta0*B0**2*n_) + 0.01 # need +0.01 for 1/T in PVMEC

                eps_a = r0/R0
                mach = np.sqrt(self['rho0']/self['rho_avg'])*self['omega0']/np.sqrt((beta0/eps_a**2)) # assumes Omega is normalised as Omega/(w_A_bar*eps_a)
                # mach = self['omega0']/np.sqrt((beta0/eps_a**2)) # assumes Omega is normalised as Omega/(w_A0*eps_a)
                self.params['mach'] = mach

            elif self['profile'] == 'rot_0': # this is the 2013 PPCF case with rot1=0
                Omega = .5*(1 + np.tanh((rstep**2 - s2)/drstep**2))
                nu_n = 2 # for comparison with reaaaally old sims
                n_ = n0*(1.-s**nu_n)
                T = np.ones_like(s)
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
            print('=== POLYFITTING FLOW ===')
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
            print('=== POLYFITTING PRESSURE ===')
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
            print('=== POLYFITTING TEMPERATURE ===')
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
        qr = self['rationalm']/(-self['n']) # n is given as a negative input because of Fourier expansion convention
        rs = self['rs'] # set to 0 to get qs = 1, this is r where q = qr, equivalent to r1 in Tom's IK work
        q0 = qr + self['delq']
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
            print('=== POLYFITTING Q-PROFILE ===')
            if C.Grid.LRFP == 'F':
                AI = np.polyfit(s2,-1./q,11)[::-1]
            elif C.Grid.LRFP == 'T':
                AI = np.polyfit(s2,-q,11)[::-1]
            else:
                print ('Insert a valid value for LRFP')
                exit()
            C.Current.AI = AI

        if self['mode_type'] == 'KH': 
            if rs == 0:
                qstep = q0+qs*self['rstep']**nu_q
            else:
                qstep = 1 - (1-q0)*(1 - (self['rstep']/rs)**nu_q)
            self.qstep = qstep

        self.dico_vmec = {'ah':AH,'ah_aux_s':AH_AUX_S,'ah_aux_f':AH_AUX_F,'at':AT,'at_aux_s':AT_AUX_S,'at_aux_f':AT_AUX_F}
        
        #Change some control parameters
        #======================================================================
        # C.Control.PREC2D_THRESHOLD = 1.0E-13 
        # C.Control.NITER_ARRAY = [1999, 3999, 3999, 3999, 8999, 8999, 8999, 8999, 8999, 25999, 39999, 99999, 129999]
        # C.Control.NS_ARRAY    = [25, 73, 211, 321, 435, 449, 463, 475, 481, 483, 485, 487, 489]
        # C.Control.FTOL_ARRAY  = [1.0e-09, 1.0e-09, 5.0e-10, 5.0e-10, 5.0e-10, 1.0e-10, 5.0e-11, 5.0e-11, 5.0e-11, 5.0e-11, 1.0e-11, 1.0e-11, 5.0E-12]
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
        print("===== RUNNING VENUS... =====")
        
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
        RationalM = self['rationalm']
        Sidebands = self['sidebands']
        stab.grid.Mmin = RationalM-Sidebands
        stab.grid.Mmax = RationalM+Sidebands
        stab.grid.Ntheta = eq.R.shape[0]
    
        stab.grid.N = 100
        stab.grid.N = eq.R.shape[1] - 2
        stab.grid.bunching = True
        stab.grid.bunchingQValues = [self.qstep - 0.01,self.qstep, self.qstep + 0.01]
        stab.grid.bunchingAmplitudes = [5.,5.,5.]
        stab.grid.bunchingSigma = [0.02,0.02,0.02]
        # stab.grid.bunchingQValues = [self.qstep]
        # stab.grid.bunchingAmplitudes = [5.]
        # stab.grid.bunchingSigma = [0.02]
        
    	#Build grid. If bunching with q values, then equilibrium quantities (radial grid s and safety factor) are required.
        stab.grid.BuildGrid(eq.s,eq.q)
    	#----------------------------------------------------------------------
        
    	#Normalize and build the equilibrium quantities in the new grid.
    	#----------------------------------------------------------------------
        eq.ChangeGrid(stab.grid.S)
        eq.Normalise()
        eq.BuildInGrid(stab.grid) # generates the derivatives of the profiles too

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
        # Omegahat0 = np.sqrt(eq.M02*2*eq.mu0*eq.P0)/(eq.B0*eps_a) # normalising to wA0
        Omegahat0 = np.sqrt(eq.M02*betahat*(np.mean(eq.rho))) # normalising to wA_bar, eq.rho is rho/rho0
        self.params['omegahat'] = Omegahat0
        # Print equilibrium quantities
        print ('Parameters at the magnetic axis:')
        print ('   Om_avg   = %.5f'%(self['omega_avg']) + ', Om_step   = %.5f'%(self['omega_step']))
        print ('   Om0hat   = %.5f'%(self['omegahat']))
        print ('   rho_avg   = %.5f'%(self['rho_avg']) + ', rho_step   = %.5f'%(self['rho_step']))
        print ('   beta_avg   = %.5f'%(self['beta_avg']) + ', beta_step   = %.5f'%(self['beta_step']))
        print ('   M0    = %.5f'%(np.sqrt(eq.M02)))
        print ('   v0/vA = %.5f'%(V0_Va))
        print ('   B0    = %.5f / %.5f [T]'%(eq.B0, self['b0']))
        print ('   R0    = %.5f / %.5f [m]'%(eq.R0, self['rmaj']))
        print ('   P0    = %.5f / %.5f [Pa]'%(eq.P0, (self['beta0']*self['b0']**2/(2.*eq.mu0)))) #P = beta0*B0**2.*n_*T/(2.*mu0*n0)
        print ('   beta0 = %.5f / %.5f %%'%(2.*eq.mu0*eq.P0/eq.B0**2., self['beta0']))
        print ('   betahat   = %.5f'%(betahat))
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
                #idx_rstep = find_nearest(stab.grid.S, self['rstep'])
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
            print ('Most unstable eigenvalue')
            print ('(Gamma/OmegaA) = %.5E + i(%.5E)'%(EV.real,EV.imag))
            
            ### Calculate peakedness
            def xi_anal(r):
                #xi_r0 = 1
                r0 = self['rstep']
                m = RationalM
                
                lam = r0**(2*m)
                
                xi = []
                for x in r:
                    if x < r0:
                        xi.append((x/r0)**(m - 1))
                    elif x >= r0:
                        xi.append((lam/(lam - 1))*((x/r0)**(m-1) - (1/lam)*(r0/x)**(m+1)))
                return xi

            # Extract the initparam
            initparam = self.scanorder[0]
            if initparam.endswith('_step'):
                initparam = initparam.replace('_step','')
            elif initparam.endswith('_avg'):
                initparam = initparam.replace('_avg','')
            elif initparam.endswith('0'):
                initparam = initparam.replace('0','')
            elif initparam.endswith('1'):
                initparam = initparam.replace('1','')
            scanprofs = {'omega':(np.asarray(eq.Omega*eq.B0/(eq.M02*eq.mu0*eq.P0)**0.5)*Omegahat0, np.asarray(eq.dOmegads*eq.B0/(eq.M02*eq.mu0*eq.P0)**0.5)*Omegahat0), 'rho':(eq.rho*self['rho0'], eq.drhods*self['rho0']), 'beta':(2*eq.P/eps_a**2, 2*eq.dPds/eps_a**2)} # dictionary to convert from inputfile names to those used in SATIRE2SFL
            # prof[0] correctly matches omegahat0 and rho0
            # Get profiles
            s = eq.s
            prof, dprofds = scanprofs[f'{initparam}']
            # Load in eigenvectors I guess
            with h5py.File(f'{self.scan_saveloc}/{self.runid}_0.h5', 'r') as f:
                xi = f['Variables']['EvaluatedModes'][f'var0_m={RationalM}'][()]
            p, p2 = self._peakedness(s, y=prof, dyds=dprofds, xi=xi, figidx=idx) # always uses the first eigenfunction for consistency i guess
            p_anal, _ = self._peakedness(s, y=prof, dyds=dprofds, xi=xi_anal(s), figidx=1)
            print(f"y0:{prof[0]}")

        output = {'EV':EV, 'v0_va':V0_Va, 'EVguess':EVguess, 'EF_file':f'{self.scan_saveloc}/{self.runid}_{idx}.h5', 'peakedness':p, 'peakedness_avgnorm':p2, 'peakedness_anal':p_anal, 'params':self.params.copy(), 'profile':self['profile'], 'shaf_diff0':shaf_diff0}
        
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
                if self['ev_guess_type'] == 'last_ev': # need a constraint so that ev_guess<0.9 always because it seems to mess up when ev_guess=1
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
        missing_scans = []
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
                    missing_scans.append(scan)
                    scans.remove(scan) # remove the skipped scan from the scanlist for the final output file
                    continue
                print(scan)
        else:
            runfile = sorted(run_saveloc.glob('*.npz'))[0]
            raw_scan = np.load(runfile, allow_pickle = True)
            rundata[runid] = raw_scan['data'].item()
            
        runinfo['scanparams'] = self.scanparams
        fixed_params = ['profile','drstep','mach','delq','rho_step','rho_avg','beta_step','beta_avg','omega_step','omega_avg','rationalm','n', 'rmaj']
        runinfo['fixedparams'] = {}
        for key in fixed_params:
            if key not in self.scanparams:
                runinfo['fixedparams'][key] = self[key]
        runinfo['scanorder'] = self.scanorder
        runinfo['scans'] = self.scans
        runinfo['missing_scans'] = missing_scans
        runinfo['timestamp'] = datetime.now().strftime("%d-%m-%y_%H:%M")
        runinfo['runid'] = self.runid
        
        fOut = f"{run_saveloc}/{self.runid}.npz"
        
        np.savez(fOut, data = rundata, info = runinfo)
        # # add any new info to the final inputfile
        # inputfile = Path.cwd() / self.inputpath
        # inputfile_nml = f90.read(inputfile)
        # inputfile_nml['info']['runid'] = self.runid
        # inputfile_nml['info']['datetime'] = datetime.now().strftime("%d-%m-%y_%H:%M")
        # with open(run_saveloc / self.inputfile, 'w') as f:
        #     f90.write(inputfile_nml, f, force=True)
        # #copy2(inputfile, run_saveloc)
        
        return # can save old runs if same input file (needs same scan parameters) and runid is provided using init_run

    def local_run(self):
        self.init_run()
        self.run(self.outpath / self.runid)
        self.save_run()
                    

        
    
    