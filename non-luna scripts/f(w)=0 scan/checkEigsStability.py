import numpy as np
import scipy as sp
from sys import getsizeof

class Stability:
    '''
		Object for stability calculations
    '''
    def __init__(self,Model):
        from VenusMHDpy.Grid import GRID
		
        self.Model = Model
        self.grid = GRID()
        self.tool = 'scipy'
	
		# Select model
        if Model == 'IdealMHD':
            from VenusMHDpy.Models import IdealMHDA as MatrixA
            from VenusMHDpy.Models import IdealMHDB as MatrixB
            self.grid.nu = [3,2,2]

        if Model == 'IdealMHDBussac':            
            from VenusMHDpy.Models import IdealMHDBussacA as MatrixA
            from VenusMHDpy.Models import IdealMHDBussacB as MatrixB
            self.grid.nu = [3,2,2]
			
        if Model == 'IdealMHDFlow-Bussac':            
            from VenusMHDpy.Models import IdealMHDFlowBussacA as MatrixA
            from VenusMHDpy.Models import IdealMHDFlowBussacB as MatrixB
            self.grid.nu = [3,2,2,3,2,2]

        if Model == 'CASTOR':  #NOT WORKING AT THE MOMENT          
            from VenusMHDpy.Models import CASTORA as MatrixA
            from VenusMHDpy.Models import CASTORB as MatrixB
            self.grid.nu = [2,3,2,2,2,2,3,3]

        if Model == 'IdealMHDEuler':            
            from VenusMHDpy.Models import IdealMHDEulerA as MatrixA
            from VenusMHDpy.Models import IdealMHDEulerB as MatrixB
            self.grid.nu = [3,2,2,2,3,3,2,2]
			
        if Model == 'IdealMHDFlow-Euler':         
            from VenusMHDpy.Models import IdealMHDFlowEulerA as MatrixA
            from VenusMHDpy.Models import IdealMHDFlowEulerB as MatrixB
            self.grid.nu = [3,2,2,2,3,3,2,2]

        if Model == 'IdealMHDperp':
            from VenusMHDpy.Models import IdealMHDperpA as MatrixA
            from VenusMHDpy.Models import IdealMHDperpB as MatrixB
            self.grid.nu = [3,2,]
			
        if Model == 'IdealMHDcyl':
            from VenusMHDpy.Models import IdealMHD_cylA as MatrixA
            from VenusMHDpy.Models import IdealMHD_cylB as MatrixB
            self.grid.nu = [2,1,1]
			
        if Model == 'LaplaceDisk':
            from VenusMHDpy.Models import LaplacianDiskA as MatrixA
            from VenusMHDpy.Models import LaplacianDiskB as MatrixB
            self.grid.nu = [3]
			
        if Model == 'SturmLiouville':
            from VenusMHDpy.Models import ModelSturmLiouvilleA_2D as MatrixA
            from VenusMHDpy.Models import ModelSturmLiouvilleB_2D as MatrixB
            self.grid.nu = [3,2]

        if Model == 'Infinite-n Ballooning':
            from VenusMHDpy.Models import Infinite_n_BallooningA_mod as MatrixA
            from VenusMHDpy.Models import Infinite_n_BallooningB_mod as MatrixB
            self.grid.nu = [3,2]
		
        self.MatrixA = MatrixA
        self.MatrixB = MatrixB
		

    def Discretize(self, eq, ntor, writeh5=False):
        from VenusMHDpy import Discretization
		
        self.ntor = ntor
		
        print ('----------------------Matrix discretisation on course------------------------')
        print ('Discretisation to be performed using %s model \n'%(self.Model))
        print ('ntor = %i'%(ntor))
        print ('mpol = [%i,%i]'%(self.grid.Mmin,self.grid.Mmax))
        print ('Ns   = %i'%(self.grid.N+2))
        print ('mtot = %i'%(self.grid.Mtot))
        if self.tool == 'scipy': self.A, self.B = Discretization.Discretize_Operators_scipy(self.grid, eq, ntor, self.MatrixA, self.MatrixB, Write=writeh5)
        print ('----------------------------------------------------------------------------- \n')
        
        
    def Solve(self,EV_guess=None,N_EV=1,maxiter=1000, which='LM', EValuesFile=None, EValuesID=0., EVectorsFile=None):
        '''
			Solve the discretised problem using scipy.sparse.linalg.eigs function.
			Input:	-A,B			: Discretised matrices of the general EV problem A*Xi = w*B*Xi.
					-EV_guess		: Eigenvalue guess for the solver.
					-N_EV			: Number of eigenvalues to calculate.
					-maxiter		: Maximum number of iterations.
					-which			: If EV_guess not specified, which eigenvalues does the solver will look for (largest magnitude, smallest real, etc. Look at scipy documentation).
					-EValuesFile	: Name of the text file with data from the simulation, incliding the eigenvalues. If None, no file is produced.
					-EValuesID		: ID for the simulation, useful for parameter scans.
					-EVectorsFile	: h5 file containing the eigenvectors. Grid is necessary to produce the file. If None, no file is produced.
        '''
        from VenusMHDpy import Solve
		
        if self.tool == 'scipy': self.Solution = Solve.SolveScipy(self.A, self.B, EV_guess, N_EV, maxiter, which, EValuesFile, EValuesID, EVectorsFile, self.grid, self.Model, self.ntor)
    
    
    def CheckZeros(self, w):
        # w is a single eigenvalue
        
        A = self.A
        B = self.B
        eigEq = A - w*B
        #print(A.shape)

        (sign, logdet) = np.linalg.slogdet(eigEq.todense())
        
        """
        try:
            det = np.log(det)
        except:
            print("det(A-omega*B) = 0, log error")
        """
        return logdet
    
