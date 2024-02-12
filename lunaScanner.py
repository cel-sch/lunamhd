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

from VenusMHDpy import SATIRE2SFL
from VenusMHDpy import Equilibriumh5
from VenusMHDpy import Stability
from VenusMHDpy import VMECInput
from VenusMHDpy.library import find_nearest

class lunaScan(object):
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
    def _readinput(self):
        # Reads input file.
        # Scanorder: list of scan parameters with order retained (it matters which is 1st)
        # Scanparams: dict of scan parameters
        self.scanparams = {}
        self.scanorder = []

        inputs = f90.read(self.inputpath)
        inputs = inputs['config_nml']

        self.scanorder = inputs['scanparams']
        if type(self.scanorder) is not list:
            self.scanorder = [self.scanorder]
        self.scandim = len(self.scanorder)

        for key in inputs:
            if key not in ['scantype','scanparams','runid']:
                if key in self.scanorder:
                    rangespecs = inputs[key]
                    paramrange = np.linspace(rangespecs[0], rangespecs[1], rangespecs[2])
                    self.scanparams[key] = paramrange
                elif key in self.params:
                    self.params[key] = inputs[key]
                elif key in self.initialisers:
                    self.initialisers[key] = inputs[key]
                else:
                    print(f"ERROR: {key} is invalid input, please review.")

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
            scandir = scandir / f'{param}_{scan[param]:.2f}'     
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
            default_inputs = f90.read(self.inputpath)
            for paramkey, paramval in scan.items(): # modify scanparams into their single values
                default_inputs['config_nml'][paramkey] = paramval
            default_inputs['config_nml']['scanparams'] = self.scanorder[0] # retain only the lowest level scan
            default_inputs.write(scan_subdir / self.inputfile, force=True) 
            
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
        
        self.label_FIXED = f'{self.runid}_{idx}' #Name for 1 point in scan e.g. test_1. VMEClabel is runid, produced when init_run is run
        
        #Read the default input file
        #C = VMECInput.ReadInputVMEC('/users/cs2427/lunamhd/VMEC/input/input.Default') # BUGFIX: will need to check if this works
        C = VMECInput.ReadInputVMEC('/home/csch/VENUS-linux/lunamhd/VMEC/input/input.Default') 
        
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
        R0 = self['rmaj']
        B0  = self['b0']
        r0 = self['a']
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
        s2 = np.linspace(0.,1.,99)
        s  = np.sqrt(s2)
        #==========================
        
        
        ### PARAMETER PROFILES ###
        #######################################################################
        
        ### Density: not a VMEC input but useful so we can set pressure wrt density.
        ### PT & rot: output AT and AH will be same as input but for the 
        # computation VMEC will normalize the profiles as T --> T/T0, 
        # Omega --> Omega/Omega0
        ### Pressure: care with P vs PVMEC

        rstep = self['rstep']
        drstep = self['drstep']

        if self['profile'] != 'rot':
            ### rot
            mach = self['mach']
            Omega = np.ones_like(s)
            Omega = (1.-s**6)
            Omega = Omega/Omega[0]
            AH = np.polyfit(s2,Omega,11)[::-1]
            # Set flow parameters to zero
            #AH = np.zeros_like(s2) # doesn't work, VMEC won't run
            
            C.Flow.AH = AH # SET FLOW PROFILE
            C.Flow.bcrit = mach # SET FLOW MAGNITUDE
            
            if self['profile'] in ['rho_t', 'rho_p']:
                ### DENSITY
                n0 = self['n0']
                n_ = .5*n0*(1 + np.tanh((rstep**2 - s2)/drstep**2))
                
                if self['profile'] == 'rho_t':
                    ### PT
                    n_ = .5*n0*(1 + np.tanh((rstep**2 - s2)/drstep**2)) + 0.05
                    
                    #T = 1/n_
                    #T = T/T[0]
                    T = .5*(1 - np.tanh((rstep**2 - s2)/drstep**2)) + 0.05
                    T = T/T[0]
                    AT = np.polyfit(s2,T,11)[::-1] # SET PT PROFILE
                    #AT = np.zeros_like(s2) # doesn't work, VMEC won't run
                    C.Flow.AT = AT
                    
                    ### PRESSURE
                    beta0 = self['beta0']
                    P = beta0*B0**2/(2*mu0*n0*T[0]) # should be constant
                    #P = beta0*B0**2/(2*mu0*n0*T[0])*(1-s2**3)
                    
                    PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
                    AM = np.polyfit(s2,PVMEC,11)[::-1]
                    C.Pressure.AM = AM # SET PRESSURE PROFILE
                    C.Pressure.PRES_SCALE = 1.
                    
                elif self['profile'] == 'rho_p':
                    ### PT
                    T = np.ones_like(s)
                    T = T/T[0]
                    AT = np.polyfit(s2,T,11)[::-1] # SET PT PROFILE
                    C.Flow.AT = AT
                    
                    ### PRESSURE
                    beta0 = self['beta0']
                    P = beta0*B0**2*n_*T/(2*mu0*n0) # stepped like n_
                    
                    PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
                    AM = np.polyfit(s2,PVMEC,11)[::-1]
                    C.Pressure.AM = AM # SET PRESSURE PROFILE
                    C.Pressure.PRES_SCALE = 1.
                    
                    # C.Pressure.PMASS_TYPE = "'cubic_spline'"
                    # C.Pressure.AM_AUX_S = s
                    # C.Pressure.AM_AUX_F = PVMEC
                    
            elif self['profile'] == 'PT':
                ### DENSITY
                n0 = self['n0']
                n_ = n0*np.ones_like(s)
                
                ### PT
                T = .5*(1 + np.tanh((rstep**2 - s2)/drstep**2)) + 0.05
                T = T/T[0]
                AT = np.polyfit(s2,T,11)[::-1] # SET PT PROFILE
                C.Flow.AT = AT
                
                ### PRESSURE
                beta0 = self['beta0']
                P0 = beta0*B0**2*n_*T/(2*mu0*n0)
                P = .5*P0*(1 + np.tanh((rstep**2 - s2)/drstep**2))
                
                PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
                AM = np.polyfit(s2,PVMEC,11)[::-1]
                C.Pressure.AM = AM # SET PRESSURE PROFILE
                C.Pressure.PRES_SCALE = 1.
                
                # C.Pressure.PMASS_TYPE = "'cubic_spline'"
                # C.Pressure.AM_AUX_S = s
                # C.Pressure.AM_AUX_F = PVMEC
            
        elif self['profile'] == 'rot':
            ### ROTATION
            mach = self['mach']
            Omega = .5*(1 + np.tanh((rstep**2 - s2)/drstep**2))
            Omega = Omega/Omega[0]
            AH = np.polyfit(s2,Omega,11)[::-1]
            
            C.Flow.AH = AH # SET FLOW PROFILE
            C.Flow.bcrit = mach # SET FLOW MAGNITUDE
            
            ### DENSITY
            n0 = self['n0']
            nu_n = self['nu_n']
            n_ = n0*(1.-s**nu_n)
            
            ### PT
            #T = np.ones_like(s)
            T = n0*(1.-s**6) + 0.01
            T = T/T[0]
            AT = np.polyfit(s2,T,11)[::-1] # SET PT PROFILE
            C.Flow.AT = AT
            
            ### PRESSURE
            beta0 = self['beta0']
            P = beta0*B0**2*n_*T/(2*mu0*n0)
            
            PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
            AM = np.polyfit(s2,PVMEC,11)[::-1]
            C.Pressure.AM = AM # SET PRESSURE PROFILE
            C.Pressure.PRES_SCALE = 1.

        self.Omega = Omega
        #======================================================================
        
        C.Current.NCURR  = 0     #0 for rotal transform, 1 for toroidal current density
        # q-profile
        #======================================================================
        qr = self['qr']
        rs = self['rs'] # set to 0 to get qs = 1, this is r where q = qr
        q0 = self['q0']
        nu_q = self['nu_q']
        
        if rs == 0:
            qs = self['qs']
        else:
            qs = (qr-q0)/rs**(nu_q)
            
        q = q0+qs*s**nu_q
        self.q = q
        
        if C.Grid.LRFP == 'F':
            AI = np.polyfit(s2,-1./q,11)[::-1]
        elif C.Grid.LRFP == 'T':
            AI = np.polyfit(s2,-q,11)[::-1]
        else:
            print ('Insert a valid value for LRFP')
            exit()
        
        C.Current.AI = AI
        
        # C.Current.PIOTA_TYPE = "'cubic_spline'"
        # C.Current.AI_AUX_S = s
        # C.Current.AI_AUX_F = q 
        
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
        DIR_VMEC = '/home/csch/VENUS-linux/lunamhd/VMEC/' # BUGFIX: MAKE SURE THIS IS SET CORRECTLY
        #DIR_VMEC = '/users/cs2427/lunamhd/VMEC/' # BUGFIX: MAKE SURE THIS IS SET CORRECTLY
        Fout = 'input.'+self.label_FIXED 
        if self['run_vmec']:
        	#Write the input file
            C.WriteInput(Fout)
        
        	#Run VMEC
            #os.system(DIR_VMEC+'./xvmec2000_flow_netcdf '+Fout)
            os.system(DIR_VMEC+'./vmec_flow '+Fout)
            print("USING VMEC VERSION vmec_flow")

            #Create the folders if not existent
            # os.system(f'mkdir -p /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/input')
            # os.system(f'mkdir -p /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/mercier')
            # os.system(f'mkdir -p /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/jxbout')
            # os.system(f'mkdir -p /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/threed1')
            # os.system(f'mkdir -p /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/wout')
            os.system(f'mkdir -p /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/input')
            os.system(f'mkdir -p /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/mercier')
            os.system(f'mkdir -p /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/jxbout')
            os.system(f'mkdir -p /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/threed1')
            os.system(f'mkdir -p /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/wout')

            #Move the output files to their folders
            # os.system('mv input.'+self.label_FIXED+f' /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/input')
            # os.system('mv wout_'+self.label_FIXED+f'.nc /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/wout')
            # os.system('mv mercier.'+self.label_FIXED+f' /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/mercier')
            # os.system('mv jxbout_'+self.label_FIXED+f'.nc /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/jxbout')
            # os.system('mv threed1.'+self.label_FIXED+f' /users/cs2427/scratch/lunamhd-data/{self.runid}/VMEC/threed1')
            os.system('mv input.'+self.label_FIXED+f' /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/input')
            os.system('mv wout_'+self.label_FIXED+f'.nc /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/wout')
            os.system('mv mercier.'+self.label_FIXED+f' /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/mercier')
            os.system('mv jxbout_'+self.label_FIXED+f'.nc /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/jxbout')
            os.system('mv threed1.'+self.label_FIXED+f' /home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/threed1')
            os.system('rm dcon_'+self.label_FIXED+'.txt')
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
        self.label_FIXED = f'{self.runid}_{idx}'
        eq = SATIRE2SFL.SATIRE2SFL(f'/home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/wout/wout_'+self.label_FIXED+'.nc')
        eq.Writeh5('eq.'+self.label_FIXED+'.h5')
        os.system('mv eq.'+self.label_FIXED+'.h5 eqFiles')
    	
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
    	
        V0_Va = np.sqrt(eq.M02*eq.mu0*eq.P0)/eq.B0
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
                EVguess = 1E-1
                EVguess = EVguess + (1.0j)*abs(n)*eq.Omega[0]
            #elif EV_guess == 'bad': #EV_guess.real < 1.0E-07 # an attempt at correcting when the EV guesses get bad
                #idx_rstep = find_nearest(stab.grid.S, self.profParams['rstep'])
                #EV_guess = 1.0E-3 + (1.0j)*abs(n)*eq.Omega[idx_rstep]
                
            print ('EV guess: %.5E + i(%.5E)'%(EVguess.real,EVguess.imag))
            #stab.Solve(EVguess,N_EV=1,EVectorsFile=f'{self.scan_saveloc}/{self.runid}_{idx}.hdf5') # runid_idx.hdf5 is where eigenvectors are stored?
            stab.Solve(EVguess,N_EV=1)
            stab.Saveh5(FileName=f'{self.scan_saveloc}/{self.runid}_{idx}')
            print ('Solution time: %.4E with N = %i' %(time.time()-t0, stab.grid.N))
    		#------------------------------------------------------------------
    
            #EV = max(stab.Solution.vals)
            EV = max(stab.vals)
            
            # Calculate peakedness
            # s2 = np.linspace(0.,1.,100) # grid
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
    		
            if self['toplot']:
                eq.plot(stab.grid, show=False)
                #stab.Solution.PlotEigenValues()
                #stab.Solution.PlotEigenVectors(eq, PlotDerivatives=False)
                stab.PlotEigenValues()
                stab.PlotEigenVectors(eq, PlotDerivatives=False)
    		#------------------------------------------------------------------

        output = {'EV':EV, 'v0_va':V0_Va, 'EVguess':EVguess, 'EF_file':f'{self.scan_saveloc}/{self.runid}_{idx}.hdf5', 'params':self.params.copy(), 'profile':self['profile']}
        
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
                self._buildVMEC(idx = vidx) # sets self.label_FIXED inside of this function
                eq = SATIRE2SFL.SATIRE2SFL(f'/home/csch/VENUS-linux/lunamhd/Output/KH/{self.runid}/VMEC/wout/wout_'+self.label_FIXED+'.nc')

            if self['run_venus']:
                # Set EV guess and calculate the growth rate
                if self['ev_guess_type'] == 'last_ev':
                    if vidx <= 2:
                        output1d[f'{scanid}_{vidx}'] = self._runVENUS(EVguess = None, idx = vidx)
                    else:
                        lastrun = output1d[f'{scanid}_{vidx-1}']
                        EVguess = lastrun['EV']
                        output1d[f'{scanid}_{vidx}'] = self._runVENUS(EVguess = EVguess, idx = vidx)
                elif self['ev_guess_type'] == 'polynom_EV':
                    ws = [] # BUGFIX: this is probably gonna be a problem for re-running scans halfway through, going to need to find a way to read results as they're being made
                    if vidx <= 3:
                        output1d[f'{scanid}_{vidx}'] = self._runVENUS(EVguess = None, idx = vidx)
                        ws.append(output1d[f'{scanid}_{vidx}']['EV'])
                    else:
                        polycoeff = vidx - 3
                        if polycoeff > 5:
                            polycoeff = 5

                        guessReal = np.polyfit(np.asarray(self.scanparams[scanparam][:vidx]),np.asarray([i.real for i in ws[:vidx]]),polycoeff)
                        guessImag = np.polyfit(np.asarray(self.scanparams[scanparam][:vidx]),np.asarray([i.imag for i in ws[:vidx]]),polycoeff)
                        
                        EVguess = np.polyval(guessReal,val) + 1j*np.polyval(guessImag,val)
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
        for n, scan in enumerate(self.scans):
            print(scan)
            scandir = self._build_scan_dir(scan)
            scanfile = sorted(scandir.glob('*.npz'))[0] # picks out the first .npz file
            scanid = f'{runid}_{n}'
            raw_scan = np.load(scanfile, allow_pickle = True)
            rundata[scanid] = raw_scan['data'].item()
        runinfo['scanparams'] = self.scanparams
        fixed_params = ['profile','drstep','mach','q0','beta0','rationalm','n', 'rmaj']
        runinfo['fixedparams'] = {}
        for key in fixed_params:
            if key not in self.scanparams:
                runinfo['fixedparams'][key] = self[key]
        runinfo['scanorder'] = self.scanorder
        runinfo['scans'] = self.scans
        runinfo['timestamp'] = datetime.now().strftime("%d-%m-%y_%H:%M")
        
        fOut = f"{run_saveloc}/{self.runid}.npz"
        
        np.savez(fOut, data = rundata, info = runinfo)
        inputpath = Path.cwd() / self.inputpath
        copy2(inputpath, run_saveloc)
        
        return # can save old runs if same input file (needs same scan parameters) and runid is provided using init_run
                    

        
    
    