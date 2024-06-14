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
        self.outputpath_root = Path('/users/cs2427/scratch/lunamhd-data')
        #self.outputpath_root = Path('/home/csch/VENUS-linux/lunamhd/Output')
        if filePath is None:
                filePath = Path(self.outputpath_root / 'KH' / f'{self.filename}')
        else:
            filePath = Path(filePath)
        self.dataFile = filePath / f'{filename}.npz'
        self._get_run_data()

    def __call__(self, variable, paramSpecs, _returnFirst = True):
        #print(f'paramSpecs={paramSpecs}, variable={variable}') # these print correctly
        # _returnFirst is a flag for whether you want call to return the first match for a variable it encounters
        # if _returnFirst = False, a list of values is returned (not yet implemented)
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

    def _make_list(self, var):
        if type(var) in [list, np.ndarray]:
            var = list(var)
        else:
            var = [var]
        return var

    def _find_nearest(self, arr, val):
        nearest_val = min(arr, key=lambda i: abs(i - val))
        return nearest_val

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
                    #if point['params'][var] != val:
                    if abs(point['params'][var] - val) > 1E-10:
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
            spar_list = self._make_list(spar_list)
            
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
        # Returns all specified variable number (varnrs) eigenfunctions for the points specified
        # varnr: the number of the eigenfunction variable being plotted (will be updated with variable names eventually)
        # scanparam: the scan parameter for which we want to see several EFs
        # spar list: if only want the variable values for certain specific scanparam values, specify those here
        # paramSpecs: any additional parameter specifications (e.g. if scan was over omega and delq, could specify to look at several delqs here) (?)  
        if spar_list is None:
            spar_list = self.info['scanparams'][scanparam] # should only load in scan parameters which are not part of fixed parameter list
        else:
            spar_list = self._make_list(spar_list)

        varnrs = self._make_list(varnrs)

        EF_dict = {}
        for p in spar_list:
            paramSpecs = deepcopy(paramSpecs)
            paramSpecs[scanparam] = p
            EF_file = self('EF_file', paramSpecs)
            EF_dict[f'{EF_file}'] = {}
            for varnr in varnrs:
                EF_dict[f'{EF_file}'][f'varnr_{varnr+1}'] = self.read_EFh5(file = EF_file, varnr = varnr)
        if _returnBoth:
            return spar_list, EF_dict
        else:
            return EF_dict

    def get_profiles_list(self, scanparam, spar_list = None, paramSpecs = {}, _returnBoth = True):
        # Returns the equilibrium profiles for the points specified
        # scanparam: the scan parameter for which we want to see several EFs
        # spar list: if only want the variable values for certain specific scanparam values, specify those here
        # paramSpecs: any additional parameter specifications (e.g. if scan was over omega and delq, could specify which delq here) (?)      
        if spar_list is None:
            spar_list = self.info['scanparams'][scanparam] # should only load in scan parameters which are not part of fixed parameter list
        else:
            spar_list = self._make_list(spar_list)

        prof_dict = {}
        for p in spar_list:
            paramSpecs = deepcopy(paramSpecs)
            paramSpecs[scanparam] = p
            EF_file = self('EF_file', paramSpecs)
            prof_dict[f'{EF_file}'] = self.read_profh5(file = EF_file)
        
        if _returnBoth:
            return spar_list, prof_dict
        else:
            return prof_dict
    
    def print_run_info(self):
        for key, val in self.info.items():
            print(f"{key}: val")

    ### PLOT FUNCTIONS ###
    def EV_plot(self, scan_specs = {}, settings = {}):
        return Plotters['Growth'](self, scan_specs=scan_specs, settings=settings)

    def EF_plot(self, varnrs, scanparam = None, spar_list = None, settings = {}):
        return Plotters['EF'](self, varnrs, scanparam, spar_list, settings)

    def profile_plot(self, scanparam = None, spar_list = None, settings = {}):
        return Plotters['Profiles'](self, scanparam, spar_list, settings)

    def multi_plot(self, readers = [], txts = {}, csvs = {}, settings = {}):
        readers.insert(0, self)
        return Plotters['Multi'](readers = readers, txts = txts, csvs = csvs, settings = settings)

    ### DATA ANALYSIS FUNCTIONS ###
    def read_EFh5(self, file, varnr = 0):
        """
        Reads h5py file for given eigenvalue run and extracts the splined eigenfunctions.
        """
        # note to self: file acts essentially as a stand-in for scanpoint
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

    def read_profh5(self, file):
        """
        Reads h5py for profile parameters.
        """
        # note to self: file acts essentially as a stand-in for scanpoint
        profile_params = {}

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

            profile_params['temp'] = 2*B0**2*T/(P0*mu0)
            profile_params['q'] = q
            profile_params['press'] = B0**2*P/mu0
            profile_params['press_rot'] = B0**2*Prot[0]/mu0
            profile_params['rho'] = rho
            if M02 == 0:
                profile_params['omega'] = np.zeros_like(s)
            else:
                profile_params['omega'] = Omega*B0/(M02*mu0*P0)**0.5
            profile_params['U'] = U/R0**2
            
            profile_params['R'] = R
            profile_params['R0'] = R0
            profile_params['Z'] = Z

            profile_params['s'] = s
        
        return profile_params




    
    