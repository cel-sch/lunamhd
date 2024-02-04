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
import csv
from VenusMHDpy import SATIRE2SFL
from VenusMHDpy import Equilibriumh5
from VenusMHDpy import Stability
from VenusMHDpy import VMECInput
from VenusMHDpy.library import find_nearest

class Solver(object):
    def __init__(self):
        self.params = {'profile':'rotation','R0':10., 'B0':1., 'r0':1., 'D':0., 'El':0., 'Tr':0., 'n':-1, 'RationalM':1,'MPOL':15, 'NTOR':0}
        # profiles = ['density_T', 'density_P', 'rotation']
        # profile: determines which parameter steps are used
        self.profParams = {'rstep':0.5, 'drstep':0.15, 'n0':1, 'nu_n':2, 'mach':0.0, 'beta0':0.05, 'qr':1., 'rs':0.3, 'q0':0.938, 'qs':None, 'nu_q':2.}
        
        self.initialisers = {'RunVMEC': True, 'RunStab': True, 'ToPlot': True, 'EVg_type': 'last_EV'}
        # EVg_type: determines what type of eigenvalue guessing system is used. 
        #'last_EV': last EV used as guess for next one, 'polynom_EV': polynomial fit used as guess for next one
        
        self.in_filepath = 'Input'
        self.in_filename = 'input.in'
        self.out_filepath = 'Output'
        self.f = self._buildOutputName() #Output filename and location
        self.VMEClabel = self.f #Label for VMEC simulation, one label used if scanning e.g. 0xffffff
    
    ### READ INPUT ###########################################################
    def _readInput(self):
        """
        Reads input file and modifies params dictionary accordingly.
        If there are parameters to scan over, adds these to scanParams dictionary.
        
        Input file only needs to contain the parameters being modified (I think).
        
        Labels:
            scan --> gets scanned over
            param --> single parameter modification
            profparam --> parameter affecting profile sizes/details
        """
        
        fInput = f'{self.in_filepath}/{self.in_filename}'
        
        self.scanParams = {}
        
        with open(fInput, 'r') as f:
            for line in f:
                data = line.split(' = ')
                if data[0] == 'scan':
                    rangeSpecs = eval(data[2])
                    paramRange = np.linspace(rangeSpecs[0], rangeSpecs[1], rangeSpecs[2])
                    self.scanParams[f'{data[1]}'] = paramRange
                elif data[0] == 'vscan':
                    scanRanges = data[2].split(' + ')
                    fullParamRange = np.array([])
                    for scanRange in scanRanges:
                        rangeSpecs = eval(scanRange)
                        paramRange = np.linspace(rangeSpecs[0], rangeSpecs[1], rangeSpecs[2])
                        fullParamRange = np.concatenate((fullParamRange, paramRange), axis=None)
                    self.scanParams[f'{data[1]}'] = fullParamRange
                if data[0] == 'profparam':
                    self.profParams[f'{data[1]}'] = eval(data[2])
                if data[0] == 'param':
                    self.params[f'{data[1]}'] = eval(data[2])
            
    ### VMEC ##################################################################                    
    def _buildVMEC(self, labelnr = 0):
        
        self.label_FIXED = f'{self.VMEClabel}_{labelnr}' #Name for 1 run in scan e.g. 0xffffff_1
        
        #Read the default input file
        C = VMECInput.ReadInputVMEC('VMEC/input/input.Default')
        
        #Modify some grid and control parameters, these get written to VMEC input
        #======================================================================
        C.Grid.MPOL = self.params['MPOL']    #Number of poloidal modes used
        C.Grid.NTOR = self.params['NTOR']    #Set 2D	
        C.Grid.LASYM = 'F'	#Weather or not to violate stellarator symmetry.
        C.Grid.NZETA = 16   #Number ot toroidal planes. Important in free boundary calculations
        C.Grid.LRFP = 'F'   #Weather to use toroidal (F) or poloidal (T) normalized flux as radial variable. Note that if poloidal flux is used, 'q' needs to be provided instead of 'iota'.
        C.FreeB.LFREEB = 'F'    #Set the simulation to be fixed boundary
        #======================================================================
        
        #Boundary
        #======================================================================
        mu0 = 4.*np.pi*1.0E-07
        e = 1.60217663E-19
        R0 = self.params['R0']
        B0  = self.params['B0']
        r0 = self.params['r0']
        D = self.params['D']  #Shafranov Shift
        El = self.params['El'] #Elongation
        Tr = self.params['Tr'] #Triangularity limit 0.08
        
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
        ### Temperature & rotation: output AT and AH will be same as input but for the 
        # computation VMEC will normalize the profiles as T --> T/T0, 
        # Omega --> Omega/Omega0
        ### Pressure: care with P vs PVMEC

        rstep = self.profParams['rstep']
        drstep = self.profParams['drstep']

        if self.params['profile'] != 'rotation':
            ### ROTATION
            mach = self.profParams['mach']
            Omega = np.ones_like(s)
            Omega = (1.-s**6)
            Omega = Omega/Omega[0]
            AH = np.polyfit(s2,Omega,11)[::-1]
            
            C.Flow.AH = AH # SET FLOW PROFILE
            C.Flow.bcrit = mach # SET FLOW MAGNITUDE
            
            if self.params['profile'] in ['density_T', 'density_P']:
                ### DENSITY
                n0 = self.profParams['n0']
                n_ = .5*n0*(1 + np.tanh((rstep**2 - s2)/drstep**2))
                
                if self.params['profile'] == 'density_T':
                    ### TEMPERATURE
                    n_ = .5*n0*(1 + np.tanh((rstep**2 - s2)/drstep**2)) + 0.05
                    
                    #T = 1/n_
                    #T = T/T[0]
                    T = .5*(1 - np.tanh((rstep**2 - s2)/drstep**2)) + 0.05
                    T = T/T[0]
                    AT = np.polyfit(s2,T,11)[::-1] # SET TEMPERATURE PROFILE
                    C.Flow.AT = AT
                    
                    ### PRESSURE
                    beta0 = self.profParams['beta0']
                    P = beta0*B0**2/(2*mu0*n0*T[0]) # should be constant
                    #P = beta0*B0**2/(2*mu0*n0*T[0])*(1-s2**3)
                    
                    PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
                    AM = np.polyfit(s2,PVMEC,11)[::-1]
                    C.Pressure.AM = AM # SET PRESSURE PROFILE
                    C.Pressure.PRES_SCALE = 1.
                    
                elif self.params['profile'] == 'density_P':
                    ### TEMPERATURE
                    T = np.ones_like(s)
                    T = T/T[0]
                    AT = np.polyfit(s2,T,11)[::-1] # SET TEMPERATURE PROFILE
                    C.Flow.AT = AT
                    
                    ### PRESSURE
                    beta0 = self.profParams['beta0']
                    P = beta0*B0**2*n_*T/(2*mu0*n0) # stepped like n_
                    
                    PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
                    AM = np.polyfit(s2,PVMEC,11)[::-1]
                    C.Pressure.AM = AM # SET PRESSURE PROFILE
                    C.Pressure.PRES_SCALE = 1.
                    
                    # trying to implement spline leads to 'Runtime Error: Factor is exactly singular' when trying to calculate EV
                    #C.Pressure.PMASS_TYPE = "'cubic_spline'"
                    #C.Pressure.AM_AUX_S = s # must match existing grid, I think
                    #C.Pressure.AM_AUX_F = PVMEC 
            
        elif self.params['profile'] == 'rotation':
            ### ROTATION
            mach = self.profParams['mach']
            Omega = .5*(1 + np.tanh((rstep**2 - s2)/drstep**2))
            Omega = Omega/Omega[0]
            AH = np.polyfit(s2,Omega,11)[::-1]
            
            C.Flow.AH = AH # SET FLOW PROFILE
            C.Flow.bcrit = mach # SET FLOW MAGNITUDE
            
            ### DENSITY
            n0 = self.profParams['n0']
            nu_n = self.profParams['nu_n']
            n_ = n0*(1.-s**nu_n)
            
            ### TEMPERATURE
            T = np.ones_like(s)
            T = T/T[0]
            AT = np.polyfit(s2,T,11)[::-1] # SET TEMPERATURE PROFILE
            C.Flow.AT = AT
            
            ### PRESSURE
            beta0 = self.profParams['beta0']
            P = beta0*B0**2*n_*T/(2*mu0*n0)
            
            PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
            AM = np.polyfit(s2,PVMEC,11)[::-1]
            C.Pressure.AM = AM # SET PRESSURE PROFILE
            C.Pressure.PRES_SCALE = 1.
        self.Omega = Omega
        #======================================================================
        
        C.Current.NCURR  = 0     #0 for rotational transform, 1 for toroidal current density
        # q-profile
        #======================================================================
        qr = self.profParams['qr']
        rs = self.profParams['rs'] # set to 0 to get qs = 1, this is r where q = qr
        q0 = self.profParams['q0']
        nu_q = self.profParams['nu_q']
        
        if rs == 0:
            qs = self.profParams['qs']
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
        DIR_VMEC = 'VMEC/'
        Fout = 'input.'+self.label_FIXED 
        if self.initialisers['RunVMEC']:
        	#Write the input file
        	C.WriteInput(Fout)
        
        	#Run VMEC
        	os.system(DIR_VMEC+'./xvmec2000_flow_netcdf '+Fout)
        	# os.system(DIR_VMEC+'./xvmec2000_netcdf '+Fout)
        
        	#Create the folders if not existent
        	os.system(f'mkdir -p VMEC/{self.VMEClabel}/input')
        	os.system(f'mkdir -p VMEC/{self.VMEClabel}/mercier')
        	os.system(f'mkdir -p VMEC/{self.VMEClabel}/jxbout')
        	os.system(f'mkdir -p VMEC/{self.VMEClabel}/threed1')
        	os.system(f'mkdir -p VMEC/{self.VMEClabel}/wout')
        
        	#Move the output files to their folders
        	os.system('mv input.'+self.label_FIXED+f' VMEC/{self.VMEClabel}/input')
        	os.system('mv wout_'+self.label_FIXED+f'.nc VMEC/{self.VMEClabel}/wout')
        	os.system('mv mercier.'+self.label_FIXED+f' VMEC/{self.VMEClabel}/mercier')
        	os.system('mv jxbout_'+self.label_FIXED+f'.nc VMEC/{self.VMEClabel}/jxbout')
        	os.system('mv threed1.'+self.label_FIXED+f' VMEC/{self.VMEClabel}/threed1')
        	os.system('rm dcon_'+self.label_FIXED+'.txt')
        #======================================================================
        
    ### VENUS #################################################################
    def _runVENUS(self, EV_guess = None, idx = 0, labelnr = 0):
        
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
        self.label_FIXED = f'{self.VMEClabel}_{labelnr}'
        eq = SATIRE2SFL.SATIRE2SFL(f'VMEC/{self.VMEClabel}/wout/wout_'+self.label_FIXED+'.nc')
        eq.Writeh5('eq.'+self.label_FIXED+'.h5')
        os.system('mv eq.'+self.label_FIXED+'.h5 eqFiles')
    	
    	#Create the stability object
        stab = Stability.Stability('IdealMHDFlow-Euler')
        eq.kappa = 0.
        
    	#Modify the default grid
    	#----------------------------------------------------------------------
        n = self.params['n'] # toroidal mode number, <0 because of how vars are expanded in n, m
        RationalM = self.params['RationalM']
        Sidebands = 5
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
        print ('   B0    = %.5f / %.5f [T]'%(eq.B0, self.params['B0']))
        print ('   R0    = %.5f / %.5f [m]'%(eq.R0, self.params['R0']))
        print ('   P0    = %.5f / %.5f [Pa]'%(eq.P0, (self.profParams['beta0']*self.params['B0']**2/(2.*eq.mu0)))) #P = beta0*B0**2.*n_*T/(2.*mu0*n0)
        print ('   beta0 = %.5f / %.5f %%'%(2.*eq.mu0*eq.P0/eq.B0**2., self.profParams['beta0']))
    	#----------------------------------------------------------------------
        
        if True:
    		# Discretize the Operators.
    		#------------------------------------------------------------------
            stab.Discretize(eq, n)
    		#------------------------------------------------------------------
    
    		# Solve
    		#------------------------------------------------------------------
            t0 = time.time()
            if EV_guess == None:
                idx_rstep = find_nearest(stab.grid.S, self.profParams['rstep'])
                #EV_guess = 1.0E-1 + (1.0j)*abs(n)*eq.Omega[idx_rstep]
                gam_guess = 1E-1
                EV_guess = gam_guess + (1.0j)*abs(n)*eq.Omega[0]
            #elif EV_guess == 'bad': #EV_guess.real < 1.0E-07
                #idx_rstep = find_nearest(stab.grid.S, self.profParams['rstep'])
                #EV_guess = 1.0E-3 + (1.0j)*abs(n)*eq.Omega[idx_rstep]
                
            print ('EV guess: %.5E + i(%.5E)'%(EV_guess.real,EV_guess.imag))
            stab.Solve(EV_guess,N_EV=1,EVectorsFile=f'{self.out_filepath}/{self.VMEClabel}_{labelnr}.hdf5')
            print ('Solution time: %.4E with N = %i' %(time.time()-t0, stab.grid.N))
    		#------------------------------------------------------------------
    
            EV = max(stab.Solution.vals)            
    		
            print ('Most unstable eigenvalue')
            print ('(Gamma/OmegaA) = %.5E + i(%.5E)'%(EV.real,EV.imag))
    		
            if self.initialisers['ToPlot']:
                eq.plot(stab.grid, show=False)
                stab.Solution.PlotEigenValues()
                stab.Solution.PlotEigenVectors(eq, PlotDerivatives=False)
    		#------------------------------------------------------------------
            
        
        return EV, V0_Va, EV_guess
    
    def doScan(self):
        
        if len(self.scanParams.keys()) == 1:
            key = self.scanParams.keys() # key being scanned over
            key = list(key)[0]
            val = self.scanParams.values() # value range of key being scanned over
            val = list(val)[0]
            
            ws = np.full(len(self.scanParams[key]), None)
            vs = np.full(len(self.scanParams[key]), None)
            pkeds = np.full(len(self.scanParams[key]), None)
            
            for idx, v in enumerate(val): #idx starts at 0
                
                # Update values of parameter being scanned over
                if key in self.profParams:
                    self.profParams[key] = v
                elif key in self.params:
                    self.params[key] = v
                
                # Run VMEC
                self._buildVMEC(labelnr = idx) # Sets label to e.g. 0xffffff and labelnr to 0, 1, 2..., produces a VMEC output file
                self.label_FIXED = f'{self.VMEClabel}_{idx}'
                eq = SATIRE2SFL.SATIRE2SFL(f'VMEC/{self.VMEClabel}/wout/wout_'+self.label_FIXED+'.nc')
                
                if self.initialisers['RunStab']:
                    # Calculate growth rate
                    if self.initialisers['EVg_type'] == 'last_EV':
                        if idx >= 2:
                            EV_guess = ws[idx-1]
                            EV_guess += 1E-3*(EV_guess.real + 1j*EV_guess.imag)
                        else:
                            EV_guess = None # set default EV_guess
                    elif self.initialisers['EVg_type'] == 'polynom_EV':
                        if idx >= 3:
                            polycoeff = idx - 3
                            if polycoeff > 5:
                                polycoeff = 5
                            
                            guessReal = np.polyfit(np.asarray(val[:idx]),np.asarray([i.real for i in ws[:idx]]),polycoeff)
                            guessImag = np.polyfit(np.asarray(val[:idx]),np.asarray([i.imag for i in ws[:idx]]),polycoeff)
                            
                            EV_guess = np.polyval(guessReal,v) + 1j*np.polyval(guessImag,v)
                            EV_guess += 1j*EV_guess.imag*1E-3
                            
                            if EV_guess.real < 1E-7 or EV_guess.real > 1E0:
                                EV_guess = ws[idx-1]  
                                EV_guess += 1E-3*(EV_guess.real + 1j*EV_guess.imag) # if the polyfit is struggling, take last calculated growth rate as guess for current one and move a little to the right
                        else:
                            EV_guess = None
                        
                    sol = self._runVENUS(EV_guess, idx, labelnr=idx)
                    w = sol[0]
                    v0va = sol[1]
                    pkedness = sol[2]
                    ws[idx] = w
                    vs[idx] = v0va
                    pkeds[idx] = pkedness
        
        else:
            print("Script currently only supports 1D scans. For 0D scans use doSweep.")
            
        return np.array([ws, vs, pkeds], dtype=object)
    
    def doSweep(self):
        if len(self.scanParams.keys()) == 0:
            
            self._buildVMEC()
            
            ### EV Guess
            EV_guess = None
            
            ### Get soln
            sol = self._runVENUS(EV_guess)
            w = sol[0]
            v0va = sol[1]
            pked = sol[2]
            
        else:
            print("Input file contains scan parameters. Use doScan.")
             
        return np.array([w, v0va, pked], dtype=object)
    
    ### BUILD OUTPUT ##########################################################
    def _buildOutputName(self):
        
        ### Randomly generate name for run
        color = random.randrange(0,2**24)
        out_filename = hex(color)
        
        return out_filename
    
    def _buildOutput(self):
        """
        Builds output file which details in filename what parameter values are.
        
        Stored in output file:
            - eigenvalue(s) for equilibri(a)um
        
        """
        
        # Make out_filepath if it doesn't exist
        os.system(f'mkdir -p {self.out_filepath}')
        
        fOutput = f'{self.out_filepath}/{self.f}.npz'
        out_filename = self.f
        
        if len(self.scanParams.keys()) == 0:
            
            scanLabel = '0D'
            output = self.doSweep()
            ws = output[0]
            vs = output[1]
            wsguess = output[2]
            
            scanParam = None
            scanVals = None
            
        elif len(self.scanParams.keys()) == 1:
            key = self.scanParams.keys() # key being scanned over
            key = list(key)[0]
            val = self.scanParams.values() # value range of key being scanned over
            val = list(val)[0]
            
            scanLabel = '1D'
            output = self.doScan()
            ws = output[0]
            vs = output[1]
            wsguess = output[2]
            
            scanParam = key
            scanVals = val
            
        ### Write output info to .csv file
        outputGrid = 'Output/outputs.csv'
        
        # Create new dictionary holding all the run information to be included in the outputs grid
        self.outparams = self.params
        self.outparams.update(self.profParams)
        self.outparams['scanlabel'] = scanLabel
        self.outparams['scanparams'] = scanParam
        self.outparams['scanvals'] = scanVals
        self.outparams['time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        csv_columns = list(self.outparams.keys())
        self.outparams['ID'] = out_filename
        self.outparams['filepath'] = self.out_filepath
        csv_columns.insert(0,'ID')
        csv_columns.insert(1, 'filepath')
        
        np.savez(fOutput, eigenvals = ws, v0vas = vs, eigenguesses = wsguess, scanvals = scanVals, scanparams = scanParam, params = self.params, profparams = self.profParams, scanlabel = scanLabel)
        
        if os.path.isfile(outputGrid): # checks whether the csv file of outputs exists already
            # check whether there are new parameters being saved in the output file --> new headers need to be written
            with open(outputGrid, 'r', newline='') as f:
                reader = csv.DictReader(f)
                headers = reader.fieldnames
                headers = list(headers)
                
            with open(outputGrid, 'a', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=csv_columns)
                headercheck = [i for i in csv_columns if i not in headers]
                headercheck += [i for i in headers if i not in csv_columns]
                #print(headercheck) prints ['profile', 'rstep', 'drstep', 'qs', 'nu_omega'], not sure why
                if len(headercheck) > 0:
                    writer.writeheader()
                writer.writerow(self.params)
        else:
            with open(outputGrid, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=csv_columns)
                writer.writeheader()
                writer.writerow(self.params)
                
        return fOutput

    ### READ OUTPUT ###########################################################                 
    def getData(self, dataFile = None, runSol = True):
        """
        Runs the solvers and unpacks the output file.
        
        Parameters:
            runSol - Determines whether or not a scan is performed when calling 
            getData function
            dataFile - If not running the solver, takes a file as input to load 
            the data. dataFile should include file location.
        """
        
        if dataFile is None:
            self._readInput()        
            data = self._buildOutput()
        else:
            data = dataFile
            
        data = np.load(data, allow_pickle = True)
        
        return data
        
        
    
    