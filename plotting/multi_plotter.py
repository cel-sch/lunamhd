# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap
from numpy import sqrt, loadtxt, float64
from pathlib import Path

from matplotlib.pyplot import subplots, show, ion, axes, tight_layout
from matplotlib.widgets import Slider, Button

import AEmhd
import lunamhd

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenval', # ['eigenval','margin_stab']
                    'x_axis_type':'initparam',
                    'rot_axis_type':'mach0', # ['mach0', 'mach1', 'omega', 'Omega', 'omegahat']
                    'axis_labels':{},
                    'reader_labels':[],
                    'AE_visible':{'gam':True, 'a_gam':False},
                    'fig_type':'general', # ['paper', 'general', 'paper_notitle']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12},
                                 'paper_notitle':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[8.5,6],
                                'paper':[4.5,3.5],
                                'paper_notitle':[4.5,3]},
                    'linestyles':{'plain':'-','asy':'D-','txt':'v--'},
                    'markersize':2,
                    'visible':{'suptitle':True, 'legend':True, 'grid':True}}

class plot_multi(object):
    def __init__(self, readers = [], txts = {}, csvs = {}, scan_specs = {}, settings = {}):
        self.readers = readers
        self.txts = txts
        self.csvs = csvs
        self.settings = {}
        defaults = deepcopy(default_settings)
        
        for key in defaults:
                self.settings[key] = defaults[key]        
        for key in settings:
            if key not in defaults:
                print(f'ERROR: {key} not found.')
            else:
                self.settings[key] = settings[key]
        
        self.labels = {'omega':'$\hat{Ω}$','beta':'$β_0$','delq':'$\hat{Δq}$','rho':'ρ','sig1':'σ', 'mach':'$\mathcal{M}$',
                        'omega0':'$\hat{Ω}_0$','omega1':'$\hat{Ω}_1$','omega_avg':'$\hat{Ω}_{avg}$', 'omega_step':'$\hat{Ω}_{step}$',
                        'rho0':'$ρ_0$','rho1':'$ρ_1$','rho_avg':'$ρ_{avg}$','rho_step':'$ρ_{step}$',
                        'beta0':'$\hat{β}_0$','beta1':'$\hat{β}_1$','beta_avg':'$\hat{β}_{avg}$','beta_step':'$\hat{β}_{step}$',
                        'T_ratio':'$T_0/T_1$',
                        'rmaj':'$R_0$','b0':'$B_0$','drstep':'$Δr_{step}$', 'q0':'$q_0$', 'rationalm':'m',
                        'eps_a':'$ε_a$','Gamma':'Γ','gam':'$\hat{γ}$','asygam':'$\hat{γ}_{asy}$',
                        'y_step':'$(y_0-y_1)/2$', 'y_avg':'$y_{avg}$','y0':'$y_0$', 'y1':'$y_1$',
                        'EV':'$\hat{ω}','a_EV':'$\hat{ω}_{asy}$','wr':'$\hat{w}_r$','asywr':'$\hat{w}_r _{asy}$'}
                        # beta is for loading from AEs, is unnormalized here (eps_a normalizations are all removed for comparison to VENUS)

        #self.outpath = Path('/users/cs2427/scratch/lunamhd-data/') # for running on viking
        self.outpath = Path(f'/home/csch/VENUS-linux/lunamhd/Output/') # for running locally
                
        self.xkeys = {}
        self.ykeys = {}
        self.scan_specs = scan_specs


        self.scankeys = {}
        self.spar_lists = {}
        # self.txtdata = {} # for more complex txt file loading

        self.open_plot()
                
    def __getitem__(self, key):
        if key in self.settings:
            return self.settings[key]
        else:
            print(f"ERROR: {key} not found")
            
    def _getlabel(self, reader, varkey):
        def _get_drive_label(format, varkey):
            if format in ['y0', 'y_0']:
                return f'${varkey}_0$'
            elif format in ['y1', 'y_1']:
                return f'${varkey}_1$'
            elif format in ['ystep', 'y_step']:
                return f'$({varkey}_0-{varkey}_1)/2$'
            elif format in ['yavg', 'y_avg']:
                return f'$({varkey}_0-{varkey}_1)/2$'

        if varkey in ['y0', 'y1', 'y_avg', 'y_step']:
            profile = reader.info['fixedparams']['profile']
            if profile in ['rhoT', 'rhoP']:
                _get_drive_label(format=varkey, varkey='ρ')
            elif profile in ['PT']:
                _get_drive_label(format=varkey, varkey='T') # or beta?
            elif profile in ['rot']:
                _get_drive_label(fomrat=varkey, varkey='$\hat{Ω}$')
        elif varkey in self.labels.keys():
            label = self.labels[varkey]
        else:
            label = varkey
        return label
 
    def save_plot(self, plotfilename, plotsaveloc = None):
        if plotsaveloc is None:
            plotsaveloc = self.outpath / self.readers[0].info['runid']
        self.fig.savefig(plotsaveloc / f'{plotfilename}')
        
    def open_plot(self):
        # Creates figure and axes
        self.fig, self.ax = subplots(figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        self.fig.set_tight_layout(True)
             
        self.scans = {}
        for reader in self.readers:
            self.scans[f'{reader}'] = reader.info['scans']
            if self.scan_specs:
                if len(self.scans[f'{reader}']) > 0:
                    self.scans[f'{reader}'] = self._make_scan_list(reader)
            self._load_x_axis(reader, self['x_axis_type'], self['rot_axis_type'])
            self._load_y_axis(reader, self['y_axis_type'])
        if self['suptitle'] is None:
            self.fig.suptitle(self._load_run_info(),fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        else:
            self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])

        if self['visible']['grid']:
            self.ax.grid()
        
        ion()
        show()
        self.draw_fig()

    def _make_scan_loop(self, reader):
        # scan_specs format e.g.: {'Omega':[1,2,3], 'beta':0.5}
        # Makes the list of e.g. [{'beta':1,'delq':0.1},{'beta':1,'delq':0.2}]
        scandim = len(reader.info['scanorder'])
        
        def loop(n = scandim - 1, scandim = scandim, scanvars = {}, scans = []):
            if n == 0:
                return [{}]
            else:
                scanparam = reader.info['scanorder'][scandim - n]
                # Loop over specified scanparam values if they are specified
                if scanparam in self.scan_specs.keys():
                    for val in self.scan_specs[scanparam]:    
                        scanvars[scanparam] = val
                        if n > 1:
                            loop(n = n - 1, scandim = scandim, scanvars = scanvars)
                        else:
                            scans.append(scanvars.copy())
                # If scanparam values not specified, load the reader scanparam values
                else:
                    for val in reader.info['scanparams'][scanparam]:    
                        scanvars[scanparam] = val
                        if n > 1:
                            loop(n = n - 1, scandim = scandim, scanvars = scanvars)
                        else:
                            scans.append(scanvars.copy())
            return scans
        return loop()
    
    def _make_scan_list(self, reader):
        if type(self.scan_specs) is list:
            scans = self.scan_specs
            self._clean_scan_list(scans, reader)
        elif type(self.scan_specs) is dict:
            scans = self._make_scan_loop(reader)
            self._clean_scan_list(scans, reader)
        else:
            print("ERROR: invalid scan_specs format provided, please input as list or dictionary.")
        return scans

    def _clean_scan_list(self, scans, reader):
        # Remove any scans which did not converge
        for scan in scans:
            if scan not in reader.info['scans']:
                scans.remove(scan)
        return scans
    
    def _make_shared_fixedparams(self):
        shared_fixedparams = {}
        for reader in self.readers:
            for key, val in reader.info['fixedparams'].items():
                if self.readers[0].info['fixedparams'][key] == val:
                    shared_fixedparams[key] = val
        return shared_fixedparams

    def _load_run_info(self):
        shared_fixedparams = self._make_shared_fixedparams()
        info = ''
        for key, val in shared_fixedparams.items():
            if type(val) in ['int', float, float64]:
                info += f"{self._getlabel(self.readers[0], key)} = {val:.2f} "
            else:
                info += f"{self._getlabel(self.readers[0], key)} = {val} "
        info = "\n".join(wrap(info, 60))
        return info
        
    def _load_x_axis(self, reader, axis_type, rot_axis_type):
        if axis_type not in ['initparam']:
            print("ERROR: axis_type not found, valid types ['initparam']")
            return
        if axis_type == 'initparam':
            initparam = deepcopy(reader.info['scanorder'][0])
            if initparam in ['mach', 'Omega', 'omega']:
                xkey = self._load_rot_axis(reader, rot_axis_type)
            else:
                xkey = initparam
            
            # need to add a switch which checks whether the output file has a y0 or not (only added omega0 to luna recently)
            # if xkey.endswith('_avg'): # y0 doesnt get saved for lunamhd
            #     xkey = xkey.replace('_avg','0')
            # elif xkey.endswith('_step'):
            #     xkey = xkey.replace('_step','0')
            self.xkeys[f'{reader}'] = xkey
            self.scankeys[f'{reader}'] = initparam
            self.spar_lists[f'{reader}'] = reader.info['scanparams'][initparam]

        if 'x' in self['axis_labels'].keys():
            self._x_ax_label = self['axis_labels']['x']
        else:
            self._x_ax_label = self._getlabel(reader, xkey) # for now: just take the xlabel from the last reader to read in, probably not a good system
            # potentially more right to load in the xlabel from the first reader
        
    def _load_y_axis(self, reader, axis_type):
        if axis_type not in ['eigenval']:
            print("ERROR: axis_type not found, valid types: ['eigenval']")
            return
        if axis_type == 'eigenval':
            if type(reader) == AEmhd.Reader.AEread:
                self.ykeys[f'{reader}'] = ['EV', 'a_EV']
            elif type(reader) == lunamhd.lunaReader.lunaRead:
                self.ykeys[f'{reader}'] = 'EV'
            self._y_ax_label = self._getlabel(reader, 'gam')

            if 'y' in self['axis_labels'].keys():
                self._y_ax_label = self['axis_labels']['y']
            else:
                self._y_ax_label = self._getlabel(reader, 'gam')

    def _load_rot_axis(self, reader, axis_type):
        if axis_type not in ['mach0', 'mach1', 'omega', 'Omega']:
            print("ERROR: rot_axis_type not found, valid types: ['mach0', 'mach1', 'omega', 'Omega']")
            return
        if axis_type in ['omega', 'Omega']:
            if type(reader) == lunamhd.lunaReader.lunaRead:
                rotkey = 'omegahat'
            elif type(reader) == AEmhd.Reader.AEread: 
                rotkey = 'omega' 
        elif axis_type in ['mach0']:
            if type(reader) == AEmhd.Reader.AEread: 
                rotkey = 'mach0'
            elif type(reader) == lunamhd.lunaReader.lunaRead:
                rotkey = 'mach'
        elif axis_type in ['mach1']:
            if type(reader) == AEmhd.Reader.AEread: 
                rotkey = 'mach1'
            elif type(reader) == lunamhd.lunaReader.lunaRead:
                rotkey = 'mach'
                print("NOTE: no way of retrieving mach1 for VENUS-MHD (yet?). mach0 being plotted.")
        return rotkey

    # def _load_txt(self, txtfile):
    #     x, y = loadtxt(txtfile, unpack = True)
        
    #     self.txtdata[f'{txtfname}'] = {}
    #     self.txtdata[f'{txtfname}']['xdata'] = x
    #     self.txtdata[f'{txtfname}']['ydata'] = y
    #     self.txtdata[f'{txtfname}']['label'] = txtfname
    #     self.txtdata[f'{txtfname}']['lstyle'] = '-'
    #     return x, y

    def _load_data(self, reader, scan = {}):
        if len(self['reader_labels']) > 0:
            if len(self['reader_labels']) != len(self.readers):
                print('ERROR: if providing reader labels, a label must be provided for every reader. Not enough/too many labels provided.')
            else:
                reader_idx = self.readers.index(reader)
                self.scanlabel = self['reader_labels'][reader_idx]
        elif scan:
            self.scanlabel = [reader.info['runid']]
            for key, keyval in scan.items():
                self.scanlabel.append(f'{self._getlabel(reader, key)}={keyval}') # empty for 1D scans
            self.scanlabel = ', '.join(self.scanlabel)             
        else:
            self.scanlabel = [reader.info['runid']]
        self.lstyle = self['linestyles']['plain']
        self.asy_lstyle = self['linestyles']['asy']

        if type(reader) == lunamhd.lunaReader.lunaRead:
            if self.xkeys[f'{reader}'] == self.scankeys[f'{reader}']: # there is definitely a way to make this into a function
                x_vals, y_vals = reader.get_1d_list(self.scankeys[f'{reader}'], variable=self.ykeys[f'{reader}'], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan) # need to check what happens if paramSpecs = None
            else:
                x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.xkeys[f'{reader}'], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                y_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
            gam_vals = [i.real for i in y_vals]
            a_gam_vals = None
            # Change x_vals from mach to omegahat if needed, probably broken but i am not fixing this rn
            if self.scankeys[f'{reader}'] == 'mach' and self['rot_axis_type'] in ['omega', 'Omega']: 
                _, x_vals = reader.get_1d_list(scanparam=self.scankeys[f'{reader}'], variable=self.xkeys[f'{reader}'], paramSpecs=scan)
            
        elif type(reader) == AEmhd.Reader.AEread:
            if reader.info['scantype'] == 'full':
                if self.xkeys[f'{reader}'] == self.scankeys[f'{reader}']:
                    x_vals, y_vals = reader.get_1d_list(self.scankeys[f'{reader}'], variable=self.ykeys[f'{reader}'][0], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan) # need to check what happens if paramSpecs = None
                    asy_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'][1], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                else:
                    x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.xkeys[f'{reader}'], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                    y_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'][0], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                    asy_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'][1], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                gam_vals = [i.imag for i in y_vals]
                a_gam_vals = [i.imag for i in asy_vals]
            elif reader.info['scantype'] == 'asy':
                if self.xkeys[f'{reader}'] == self.scankeys[f'{reader}']:
                    x_vals, asy_vals = reader.get_1d_list(self.scankeys[f'{reader}'], variable=self.ykeys[f'{reader}'][1], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan) # need to check what happens if paramSpecs = None
                else:
                    x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.xkeys[f'{reader}'], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                    asy_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'][1], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                a_gam_vals = [i.imag for i in asy_vals]
            # Change x_vals from omega to mach0 or mach1, pretty sure this is broken atm because of new omega options
            if self.scankeys[f'{reader}'] == 'omega' and self['rot_axis_type'] in ['mach0', 'mach1']: 
                _, x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.xkeys[f'{reader}'], paramSpecs = scan)

        # for txtfname, txtf in self.txts:
        #     self._load_txt(txtfile = txtf)

        return x_vals, gam_vals, a_gam_vals

    def plot_vals(self, reader = None, txtfname = None, csvfname = None, scan = {}):
        if reader:
            if type(reader) == lunamhd.lunaReader.lunaRead:
                x_vals, gam_vals, _ = self._load_data(reader = reader, scan = scan)
                self.ax.plot(x_vals, gam_vals, self.lstyle, label=f'{self.scanlabel}', markersize=self['markersize'])
            elif type(reader) == AEmhd.Reader.AEread:
                if self['AE_visible']['gam']:
                    x_vals, gam_vals, _ = self._load_data(reader = reader, scan = scan)
                    self.ax.plot(x_vals, gam_vals, self.lstyle, label=f'{self.scanlabel} DE', markersize=self['markersize'])
                if self['AE_visible']['a_gam']:
                    x_vals, _, a_gam_vals = self._load_data(reader = reader, scan = scan)
                    self.ax.plot(x_vals, a_gam_vals, self.asy_lstyle, label=f'{self.scanlabel} AE', markersize=self['markersize'])
        elif txtfname:
            x_vals, y_vals = loadtxt(self.txts[f'{txtfname}'], unpack = True)
            self.ax.plot(x_vals, y_vals, '--v', label = txtfname, markersize=self['markersize'])
        elif csvfname:
            x_vals, y_vals = loadtxt(self.csvs[f'{csvfname}'], delimiter = ',', unpack = True, skiprows=1)
            y_vals = [i*0.1 for i in y_vals] # need a more general way to renormalize this
            self.ax.plot(x_vals, y_vals, '--v', label = csvfname, markersize=self['markersize'])

        
    def draw_fig(self):
        for reader in self.readers:
            scans = self.scans[f'{reader}']
            if scans:
                for scan in scans:
                    self.plot_vals(reader = reader, scan = scan)
            else:
                self.plot_vals(reader = reader)

        for txtfname in self.txts.keys():
            self.plot_vals(txtfname = txtfname)

        for csvfname in self.csvs.keys():
            self.plot_vals(csvfname = csvfname)
                            
        self.ax.set_ylabel(self._y_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis'])
        self.ax.set_xlabel(self._x_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis']) 
        if self.scanlabel:
            self.ax.legend()
            self.ax.legend_.set_visible(self['visible']['legend'])
        
        self.fig.canvas.draw_idle()
