# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 16:44:05 2022

@author: celin
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import time
from VenusMHDpy import SATIRE2SFL
from VenusMHDpy import Equilibriumh5
from VenusMHDpy import Stability
from VenusMHDpy import VMECInput


'''
	Runs Fixed boundary 2D VMEC, transforms the equilibrium to straight field line coordinates, and solves the stability problem using VENUS-MHDpy. 
'''


RunVMEC = True
RunStab = True
ToPlot  = True
label_FIXED = 'IK_Ugly'        #Label for the VMEC simulation

#----------------------------------------------------------------------------------------------#
#=======================================VMEC FIXED BOUNDARY====================================#
#----------------------------------------------------------------------------------------------#

GammaNorm = []
OmegaNorm = []
V0_VaArray = []
#param = 'mach'
ParamArray = np.linspace(0.0,0.25,20)

xaxis = r'$\mathcal{M}_0$'

k = 0 # counting parameter
for x in ParamArray:
    mach = x

    #Read the default input file
    C = VMECInput.ReadInputVMEC('VMEC/input/input.Default')
    
    #Modify some grid and control parameters
    #==========================================================================
    C.Grid.MPOL = 15    #Number of poloidal modes used
    C.Grid.NTOR = 0		#Set 2D	
    C.Grid.LASYM = 'F'	#Weather or not to violate stellarator symmetry.
    C.Grid.NZETA = 16   #Number ot toroidal planes. Important in free boundary calculations
    C.Grid.LRFP = 'F'   #Weather to use toroidal (F) or poloidal (T) normalized flux as radial variable. Note that if poloidal flux is used, 'q' needs to be provided instead of 'iota'.
    C.FreeB.LFREEB = 'F'    #Set the simulation to be fixed boundary
    #==========================================================================
    
    #Boundary
    #====================================================================
    mu0 = 4.*np.pi*1.0E-07
    e = 1.60217663E-19
    R0 = 3.33
    B0  = 1.0
    r0 = 1.0
    D = 0.0  #Shafranov Shift
    El = 0.0 #Elongation
    Tr = 0.0 #Triangularity limit 0.08
    
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
    #====================================================================
    
    #Grid
    #==========================
    s2 = np.linspace(0.,1.,100)
    s  = np.sqrt(s2)
    #==========================
    
    
    #Density profile. This is not a VMEC input, but sometimes it is useful to 
    #define one so that we can set the Pressure with respect to the density
    #==============================================================================
    n0 = 2.5E+21
    nu_n = 2.
    n_ = n0*(1.-s**nu_n)
    #==============================================================================
    
    #Temperature and rotation profiles. Note that The output AT and AH will be the 
    #same as the input, but for the computation VMEC will normalize the profiles as
    #T --> T/T0, Omega --> Omega/Omega0. 
    #==============================================================================
    #Set the temperature equal to the pressure, which makes the density constant
    #AT = np.polyfit(s2,P/P0,11)[::-1]
    #Set the temperature constant, which makes the density equal to the pressure
    T = np.ones_like(s)
    #AT = [1]
    
    T = T/T[0]
    AT = np.polyfit(s2,T,11)[::-1]
    C.Flow.AT = AT
    
    #------------------------------------------------------------------------------
    #mach = 0.5
    
    nu_omega = 2.
    Omega = 1.-s**nu_omega
    AH = np.polyfit(s2,Omega,11)[::-1]
    C.Flow.AH = AH # SET FLOW PROFILE
    C.Flow.bcrit = mach # SET FLOW MAGNITUDE
    #==============================================================================
    
    #Pressure Profile
    #===============================================================
    #Parameters for the model.
    #----------------------------------------
    beta0 = 0.054
    P0 = B0**2*beta0/(2*mu0)
    P = P0*n_*T/n0
    
    #P = beta0*B0**2.*n_*T/(2.*mu0*n0) # same thing  as above i think
    
    PVMEC = P*np.exp(-0.5*mach**2*Omega**2./T)
    
    AM = np.polyfit(s2,PVMEC,11)[::-1]
    C.Pressure.AM = AM # SET PRESSURE PROFILE
    C.Pressure.PRES_SCALE = 1.
    #----------------------------------------
    
    #Set the profiles in VMEC
    C.Pressure.AM = AM
    C.Pressure.PRES_SCALE = 1.
    #===============================================================
    
    C.Current.NCURR  = 0     #0 for rotational transform, 1 for toroidal current density
    
    # #Rotational transform
    # #===============================================================
    #Arguments
    #----------------------------------------
    qr = 1
    rs = 0.3
    q0 = 0.868
    nu_q = 2.
    
    qs = (qr-q0)/rs**(nu_q)
    q = q0+qs*s**nu_q
    
    if C.Grid.LRFP == 'F':
    	AI = np.polyfit(s2,-1./q,11)[::-1]
    elif C.Grid.LRFP == 'T':
    	AI = np.polyfit(s2,-q,11)[::-1]
    else:
    	print ('Insert a valid value for LRFP')
    	exit()
    #----------------------------------------
    
    #Set the profiles in VMEC
    C.Current.AI = AI
    # #===============================================================
    
    #Change some control parameters
    #=============================================================================================================================================
    C.Control.PREC2D_THRESHOLD = 1.0E-13 
    C.Control.NITER_ARRAY = [1999, 3999, 3999, 3999, 8999, 8999, 8999, 8999, 8999, 25999, 39999, 99999, 129999]
    C.Control.NS_ARRAY    = [25, 73, 211, 321, 435, 449, 463, 475, 481, 483, 485, 487, 489]
    C.Control.FTOL_ARRAY  = [1.0e-09, 1.0e-09, 5.0e-10, 5.0e-10, 5.0e-10, 1.0e-10, 5.0e-11, 5.0e-11, 5.0e-11, 5.0e-11, 1.0e-11, 1.0e-11, 5.0E-12]
    #=============================================================================================================================================
    
    #Run VMEC Fixed boundary VMEC
    #========================================================
    DIR_VMEC = 'VMEC/'
    Fout = 'input.'+label_FIXED 
    if RunVMEC:
    	#Write the input file
    	C.WriteInput(Fout)
    
    	#Run VMEC
    	os.system(DIR_VMEC+'./xvmec2000_flow_netcdf '+Fout)
    	# os.system(DIR_VMEC+'./xvmec2000_netcdf '+Fout)
    
    	#Create the folders if not existent
    	os.system('mkdir -p VMEC/input')
    	os.system('mkdir -p VMEC/mercier')
    	os.system('mkdir -p VMEC/jxbout')
    	os.system('mkdir -p VMEC/threed1')
    	os.system('mkdir -p VMEC/wout')
    
    	#Move the output files to their folders
    	os.system('mv input.'+label_FIXED+' VMEC/input')
    	os.system('mv wout_'+label_FIXED+'.nc VMEC/wout')
    	os.system('mv mercier.'+label_FIXED+' VMEC/mercier')
    	os.system('mv jxbout_'+label_FIXED+'.nc VMEC/jxbout')
    	os.system('mv threed1.'+label_FIXED+' VMEC/threed1')
    	os.system('rm dcon_'+label_FIXED+'.txt')
    #========================================================
    
    #Run stability code
    #=======================================================================================================
    if RunStab:
    	
       	#Read equilibrium from VMEC output file and transform it into SFL
       	eq = SATIRE2SFL.SATIRE2SFL('VMEC/wout/wout_'+label_FIXED+'.nc')
       	eq.Writeh5('eq.'+label_FIXED+'.h5')
       	os.system('mv eq.'+label_FIXED+'.h5 eqFiles')
       	
       	#Create the stability object
       	stab = Stability.Stability('IdealMHDFlow-Euler')
       	eq.kappa = 0.
       
       	#Modify the default grid
       	#---------------------------
       	n = -1
       	RationalM = 1
       	Sidebands = 5
       	stab.grid.Mmin = RationalM-Sidebands
       	stab.grid.Mmax = RationalM+Sidebands
       	stab.grid.Ntheta = eq.R.shape[0]
       
       	stab.grid.N = 100
       	stab.grid.bunching = True
       	stab.grid.bunchingQValues = [1.0,1.1,1.2]
       	stab.grid.bunchingAmplitudes = [5.,5.,5.]
       	stab.grid.bunchingSigma = [0.02,0.02,0.02]
       	
       	#Build grid. If bunching with q values, then equilibrium quantities (radial grid s and safety factor) are required.
       	stab.grid.BuildGrid(eq.s,eq.q)
       	#--------------------------------------			
       	
       	
       	#Normalize and build the equilibrium quantities in the new grid.
       	#-----------------------------------------------
       	eq.ChangeGrid(stab.grid.S)
       	eq.Normalise()
       	eq.BuildInGrid(stab.grid)
       	
       	#eq.plot(stab.grid, show=True)
       	
       	V0_Va = np.sqrt(eq.M02*eq.mu0*eq.P0)/eq.B0
        V0_VaArray.append(V0_Va)
       	print ('Parameters at the magnetic axis:')
       	print ('   M0    = %.5f'%(np.sqrt(eq.M02)))
       	print ('   v0/vA = %.5f'%(V0_Va))
       	print ('   B0    = %.5f'%(eq.B0))
       	print ('   R0    = %.5f'%(eq.R0))
       	print ('   P0    = %.5f'%(eq.P0))
       	print ('   beta0 = %.5f'%(2.*eq.mu0*eq.P0/eq.B0**2.))
       	#-----------------------------------------------
        
        if True:
    		# Discretize the Operators.
    		#----------------------------------------------------------------------------
            stab.Discretize(eq, n)
      		#----------------------------------------------------------------------------
      
      		# Solve
      		#--------------------------------------------------------------------
            t0 = time.time()
            #Set eigenvalue guess based on trend
            #----------------------------------------------
            if k == 0:
                EV_guess = 1.0E-03 + (1.0j)*abs(n)*(eq.M02*eq.mu0*eq.P0/eq.B0**2.)**0.5
            else:
                EV_guess = GammaNorm[-1]
            #----------------------------------------------
            print ('EV guess: %.5E + i(%.5E)'%(EV_guess.real,EV_guess.imag))
            stab.Solve(EV_guess,N_EV=1, EValuesFile="test.txt", EValuesID=f"{x}")
            # stab.Solve(self,EV_guess=None,N_EV=1,maxiter=1000, which='LM', Values_txt=None, ValuesID=0., Vectors_h5=None, grid=None):
            print ('Solution time: %.4E with N = %i' %(time.time()-t0, stab.grid.N))
            #--------------------------------------------------------------------
            
            EV = max(stab.Solution.vals)
            #if EV < 1.0E-06:
                #break
            GammaNorm.append(EV)
            
            print ('Most unstable eigenvalue')
            print ('(Gamma/OmegaA) = %.5E + i(%.5E)'%(EV.real,EV.imag))
            k =+ 1
            		
            
            if ToPlot:
                eq.plot(stab.grid, show=False)
                #stab.Solution.PlotEigenValues()
                #stab.Solution.PlotEigenVectors(eq, PlotDerivatives=False)
            		#-------------------------------------------------------------------
    	#=======================================================================================================

# Saving data
path = "Output/IK"
filename = "IK_Ugly.txt"
datapath = f"{path}/{filename}"

data = np.column_stack([V0_VaArray, ParamArray, GammaNorm])
np.savetxt(datapath, [], header="v0/va, Mach, γ/ωA")
with open(datapath, 'a') as f:
    np.savetxt(f, data)
    
