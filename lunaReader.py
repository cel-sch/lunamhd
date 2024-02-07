import h5py
import numpy as np
import scipy.integrate as spi
from pathlib import Path
from copy import deepcopy

from VenusMHDpy.Bsplines import Bspline
from VenusMHDpy.Grid import GRID
from lunamhd.plotting import Plotters

class lunaRead(object):
    def __init__(self, filename, filePath = None):
        self.filename = filename
        if filePath is None:
                filePath = Path(f'/users/cs2427/scratch/lunamhd-data/KH/{filename}')
        else:
            filePath = Path(filePath)
        self.dataFile = filePath / f'{filename}.npz'
        self._get_run_data()

    def __call__(self, variable, paramSpecs):
        #print(f'paramSpecs={paramSpecs}, variable={variable}') # these print correctly
        scankey, pointkey = self.get_point_label(paramSpecs=paramSpecs)
        if pointkey:
            scan = self.data[scankey]
            for key in scan[pointkey].keys():
                if key == variable:
                    return scan[pointkey][key]
                elif type(scan[pointkey][key]) == dict:
                    for skey in scan[pointkey][key].keys():
                        if skey == variable:
                            return scan[pointkey][key][skey]
            print("ERROR: cannot find variable")
            return
        else:
            return

    def _get_run_data(self, dataFile = None):
        if dataFile is None:
            dataFile = self.dataFile
        else:
            self.dataFile = dataFile
        raw_data = np.load(dataFile, allow_pickle = True)
        info = raw_data['info'].item() # contains scanparams, timestamp
        data = raw_data['data'].item() # {scan1:{data:{point1, point2, ...},info:info}, ...}
        self.data = data  # data isn't building correctly
        self.info = info
        return

    def get_point_label(self, paramSpecs):
        # NOTE: this relies on ordered dictionaries and will not work for python versions older than 3.6
        # Currently only returns the first point matching the paramSpecs given.
        initparam = self.info['scanorder'][0]
        if initparam in paramSpecs.keys(): # move initparam to beginning of dict
            paramSpecs = {f'{initparam}':paramSpecs.pop(initparam), **paramSpecs}
            
        for scankey, scan in self.data.items():
            for pointkey, point in scan.items():
                isrun = True
                for var, val in paramSpecs.items():
                    if point['params'][var] != val:
                        isrun = False
                    #else: isrun = True breaks things to always pick the same point, unsure why
                if isrun:
                    return scankey, pointkey
        if not isrun:
            print("ERROR: Could not find run")
            return None

    def get_1d_list(self, scanparam, variable, spar_list = None, paramSpecs = {}, _returnBoth = True):
        # Returns a 1D list of variable values 
        # variable: the variable values of which are being looked up
        # scanparam: the value over which was scanned
        # spar list: if only want the variable values for certain specific scanparam values, specify those here
        # paramSpecs: any additional parameter specifications (e.g. if scan was over omega and delq, could specify which delq here) (?)
        if spar_list is None:
            spar_list = self.info['scanparams'][scanparam] # should only load in scan parameters which are not part of fixed parameter list
        else:
            spar_list = list(spar_list)
            
        var_list = []
        for p in spar_list:
            paramSpecs = deepcopy(paramSpecs)
            paramSpecs[scanparam] = p 
            var_list.append(self(variable, paramSpecs))
        if _returnBoth:
            return spar_list, var_list
        else:
            return var_list

    def get_eigenfunc_list(self, varnrs, scanparam, spar_list = None, paramSpecs = {}, _returnBoth = True):
        # Returns EF values
        # varnr: the number of the eigenfunction variable being plotted (will be updated with variable names eventually)
        # scanparam: the value over which was scanned
        # spar list: if only want the variable values for certain specific scanparam values, specify those here
        # paramSpecs: any additional parameter specifications (e.g. if scan was over omega and delq, could specify which delq here) (?)      
        if spar_list is None:
            spar_list = self.info['scanparams'][scanparam] # should only load in scan parameters which are not part of fixed parameter list
        else:
            spar_list = list(spar_list)
        varnrs = [varnrs]

        EF_file_list = []
        for p in spar_list:
            paramSpecs = deepcopy(paramSpecs)
            paramSpecs[scanparam] = p 
            EF_file_list.append(self('EF_file', paramSpecs))
            
        EF_list = {}
        if len(EF_file_list) == 1: # looking at one point
            for varnr in varnrs:
                EF_list[f'varnr_{varnr+1}'] = self.read_EFh5(file = EF_file_list[0], varnr = varnr)
        elif len(EF_file_list) > 1: # looking at multiple points
            if len(varnrs) > 1:
                print("ERROR: trying to retrieve multiple eigenfunctions for multiple scan points. Either specify one scan point or one eigenfunction to look at.")
                return
            for EF_file in EF_file_list:
                EF_list[f'{EF_file}'] = self.read_EFh5(file = EF_file, varnr = varnrs[0])

        return EF_list
    
    def print_run_info(self):
        for key, val in self.info.items():
            print(f"{key}: val")

    ### PLOT FUNCTIONS ###
    def basic_plot(self, settings = {}):
        return Plotters['Growth'](self, settings)

    ### DATA ANALYSIS FUNCTIONS ###
    def read_EFh5(self, file, varnr = 0):
        """
        Reads h5py file for given eigenvalue run and splines the eigenfunctions.
        """
        grid = GRID()
        with h5py.File(file, 'r') as f:
            ms = f['grid']['m'][()] # array of poloidal mode numbers
            nus = f['grid']['nu'][()] # grid.nu is determined by which model is being used, for IdealMHDFlow-Euler is [3,2,2,2,3,3,2,2] because 8 variables
            N = f['grid']['N'][()] # number of grid points?
            #S = f['grid']['S'][()]
            
            ### Update grid with new parameters
            grid.N = N
            grid.nu = nus
            grid.S = f['grid']['S'][()]
            
            grid.Mmin = min(f['grid']['m'][()])
            grid.Mmax = max(f['grid']['m'][()])
            
            grid.knots = f['grid']['knots'][()]
            grid.sk = f['grid']['sk'][()]
            
            # not all grid parameters are specified here, bspline uses: knots, sk, N, S
            
            ### Build Bspline arrays
            r = np.linspace(0.,1.,10000)
            BspCalc = []
            
            nu = nus[varnr]      
            vec_allms = f['variables'][f'var{varnr}'] # all poloidal mode numbers for variable varnr (in h5py vars go from 0 to 7)

            #Create the Bspline arrays
            #-------------------------------------------------------------------------------------------------
            if nu not in BspCalc:
                BspCalc.append(nu)	
                vars()['Bsp'+str(nu)] = np.zeros(shape=(r.size, N+1+nu))
                
                for j in range(N):
                    l = np.ones(len(r), dtype=int)*j-nu
                    vars()['Bsp'+str(nu)][:,j] = Bspline(r,l,nu,grid,der=0)
                    
            mixB = max(nus)-nu # 8 - order of variable (2 or 3 generally)
            
            mode_allms = {}
            for j in range(f['grid']['Mtot'][()]): # selects vector for every poloidal mode number
                vec = vec_allms[j]
                mode = np.sum(vars()['Bsp'+str(nu)]*vec, axis=1) # performs a bspline on the vector?
                mode_allms[f'm={ms[j]}'] = mode

            return mode_allms



    
    