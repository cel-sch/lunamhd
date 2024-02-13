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
                #filePath = Path(f'/users/cs2427/scratch/lunamhd-data/KH/{filename}')
                filePath = Path(f'/home/csch/VENUS-linux/lunamhd/Output/KH/{filename}')
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
        # scanparam: the scan parameter for which we want to see several EFs
        # spar list: if only want the variable values for certain specific scanparam values, specify those here
        # paramSpecs: any additional parameter specifications (e.g. if scan was over omega and delq, could specify which delq here) (?)      
        if spar_list is None:
            spar_list = self.info['scanparams'][scanparam] # should only load in scan parameters which are not part of fixed parameter list
        else:
            spar_list = list(spar_list)
        varnrs = list(varnrs)

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

    def EF_plot(self, varnrs, scanparam = None, spar_list = None, settings = {}):
        return Plotters['EF'](self, varnrs, scanparam, spar_list, settings)

    ### DATA ANALYSIS FUNCTIONS ###
    def read_EFh5(self, file, varnr = 0):
        """
        Reads h5py file for given eigenvalue run and splines the eigenfunctions.
        """
        grid = GRID()
        mode_allms = {}
        with h5py.File(file, 'r') as f:
            grid.S = f['Grid']['S'][()]
            r = grid.S
            for key in [x for x in f['Variables']['EvaluatedModes'].keys() if f'var{varnr}' in x]:
                mode = f['Variables']['EvaluatedModes'][key][()]
                mval = key.split("=",1)[1]
                mode_allms[f'm={mval}'] = mode
        return mode_allms



    
    