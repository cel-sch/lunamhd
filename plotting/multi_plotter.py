# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap
from numpy import sqrt, loadtxt, float64, pi, linspace
from pathlib import Path

from matplotlib.pyplot import subplots, show, ion, axes, tight_layout
from matplotlib.widgets import Slider, Button

import AEmhd
import lunamhd

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenval', # ['eigenval','margin_stab']
                    'y_axis_lims':None,
                    'x_axis_type':'xstep_norm', # ['initparam', 'peakedness', 'peakedness_anal', 'xstep_norm']
                    'x_axis_lims':None,
                    'rot_axis_type':'mach0', # ['mach0', 'mach1', 'omega', 'Omega', 'omegahat']
                    'axis_labels':{},
                    'reader_labels':[],
                    'AE_visible':{'gam':True, 'a_gam':False},
                    'fig_type':'paper', # ['paper', 'general', 'paper_notitle']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12},
                                 'paper_notitle':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[7.5,6],
                                'paper':[4.5,4],
                                'paper_notitle':[4.5,3]},
                    'linestyles':{'plain':'-','asy':'D-','txt':'v--'},
                    'own_ls':[],
                    'markfreq':[],
                    'markersize':4,
                    'plotrange':{'pstart':None, 'pstop':None},
                    'visible':{'suptitle':True, 'legend':True},
                    'legend_loc':'best'}

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

        self.labels = {'omega':'$\hat{Ω}$','beta':'$β$','delq':'$\hat{Δq}$','rho':'ρ','sig1':'σ', 'mach':'$\mathcal{M}$',
                        'omega0':'$\hat{Ω}_0$','omega1':'$\hat{Ω}_1$','omega_avg':'$\hat{Ω}_{avg}$', 'omega_step':'$\hat{Ω}_{step}$',
                        'rho0':'$ρ_0$','rho1':'$ρ_1$','rho_avg':'$ρ_{avg}$','rho_step':'$ρ_{step}$',
                        'beta0':'$\hat{β}_0$','beta1':'$\hat{β}_1$','beta_avg':'$\hat{β}_{avg}$','beta_step':'$\hat{β}_{step}$',
                        'T_ratio':'$T_0/T_1$',
                        'rho0_pkd':'$ρ$ peakedness', 'omega0_pkd':'$\hat{Ω}$ peakedness', 'beta0_pkd':'$\hat{β}$ peakedness',
                        'rho_xstep':'$(ρ_0-ρ_1)/ρ_{avg}$', 'omega_xstep':'$(\hat{Ω}_0-\hat{Ω}_1)/\hat{Ω}_{avg}$', 'beta_xstep':'$(\hat{β}_0-\hat{β}_1)/\hat{β}_{avg}$',
                        'rho_xstep':'$\\tilde{ρ}_{step}$', 'omega_xstep':'$\\tilde{Ω}_{step}$', 'beta_xstep':'$\\tilde{β}_{step}$',
                        'rmaj':'$R_0$','b0':'$B_0$','drstep':'$Δr_{step}$', 'q0':'$q_0$', 'rationalm':'m',
                        'eps_a':'$ε_a$','Gamma':'Γ','gam':'$\hat{γ}$','asygam':'asymptotic $\hat{γ}$',
                        'y_step':'$(y_0-y_1)/2$', 'y_avg':'$y_{avg}$','y0':'$y_0$', 'y1':'$y_1$',
                        'EV':'$\hat{ω}$','a_EV':'asymptotic $\hat{ω}$','wr':'$\hat{ω}_r$','asywr':'asymptotic $\hat{ω}_r$'}

        #self.outpath = Path('/users/cs2427/scratch/lunamhd-data/') # for running on viking
        self.outpath = Path(f'/home/csch/VENUS-linux/lunamhd/Output/KH') # for running locally
                
        self.initparams = {}
        self.xkeys = {}
        self.xkey0s = {}
        self.ykeys = {}
        self.scan_specs = scan_specs


        self.scankeys = {}
        self.spar_lists = {}
        # self.txtdata = {} # for more complex txt file loading

        self.wA_avgNorm = True # convert VENUS outputs to wA_avg normalisation or not

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
        self.fig.savefig(Path(plotsaveloc) / f'{plotfilename}')
        
    def open_plot(self):
        # Creates figure and axes
        self.fig, self.ax = subplots(figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        if self['x_axis_lims']:
            self.ax.set_xlim(self['x_axis_lims'][0],self['x_axis_lims'][1])
        if self['y_axis_lims']:
            self.ax.set_ylim(self['y_axis_lims'][0],self['y_axis_lims'][1])
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

        # if self['visible']['grid']:
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
                if key in self.readers[0].info['fixedparams'].keys() and self.readers[0].info['fixedparams'][key] == val:
                    shared_fixedparams[key] = val
                elif key in self.readers[0].info['fixedparams'].keys():
                    try:
                        del shared_fixedparams[key]
                    except:
                        pass
        return shared_fixedparams

    def _load_run_info(self): 
        shared_fixedparams = self._make_shared_fixedparams()
        info = ''
        for key, val in shared_fixedparams.items():
            if type(val) in ['int', float, float64]:
                info += f"{self._getlabel(self.readers[0], key)} = {val:.4f} "
            else:
                info += f"{self._getlabel(self.readers[0], key)} = {val} "
        info = "\n".join(wrap(info, 60))
        return info
        
    def _load_x_axis(self, reader, axis_type, rot_axis_type):
        if axis_type not in ['initparam', 'peakedness', 'peakedness_anal','xstep_norm']:
            print("ERROR: axis_type not found, valid types ['initparam', 'peakedness', 'peakedness_anal', 'xstep_norm']")
            return

        initparam = deepcopy(reader.info['scanorder'][0]) 
        self.initparams[f'{reader}'] = initparam
        if initparam in ['mach', 'Omega', 'omega']:
            xkey = self._load_rot_axis(reader, rot_axis_type)
        else:
            xkey = initparam
            xkey0 = initparam
        # need to add a switch which checks whether the output file has a y0 or not (only added omega0 to luna recently)
        if xkey.endswith('_avg'):
            xkey = xkey.replace('_avg','')
        elif xkey.endswith('_step'):
            xkey = xkey.replace('_step','')

        if axis_type in ['peakedness', 'peakedness_anal']:
            xkey = f'{xkey}_pkd'

        self.xkeys[f'{reader}'] = xkey
        self.xkey0s[f'{reader}'] = xkey0
        self.scankeys[f'{reader}'] = initparam
        self.spar_lists[f'{reader}'] = reader.info['scanparams'][initparam]

        if 'x' in self['axis_labels'].keys():
            self._x_ax_label = self['axis_labels']['x']
        else:
            self._x_ax_label = self._getlabel(reader, xkey) # for now: just take the xlabel from the last reader to read in, probably not a good system
            # potentially more right to load in the xlabel from the first reader
            if axis_type == 'xstep_norm':
                self._x_ax_label = self._getlabel(reader, f'{xkey}_xstep')
            elif axis_type == 'initparam':
                self._x_ax_label = self._getlabel(reader, xkey0)
        
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
            if self['x_axis_type'] in ['peakedness', 'peakedness_anal']:
                rotkey = f'{self.xkey}_pkd'
            else:
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

    def load_stepavg(self, reader, prof, scan = {}):
        mu0 = 4.*pi*1.0E-07
        steps = []
        avgs = []
        for idx in range(len(self.spar_lists[f'{reader}'])):
            for key, vals in reader.get_profiles_list(scanparam = self.initparams[f'{reader}'], spar_list = self.spar_lists[f'{reader}'][idx], paramSpecs = scan, _returnBoth = False).items():
                EF_file = key
                prof_dict = vals
            EF_file = EF_file.split('/')[-1]

            B0 = prof_dict['B0']
            P0 = prof_dict['P0']
            R0 = prof_dict['R0']
            M02 = prof_dict['M02']

            if prof in ['omega', 'Omega']:
                Omega = prof_dict['omega']
                Omega = Omega*B0/(M02*mu0*P0)

                step = (Omega[0] - Omega[-3])/2
                avg = (Omega[0] + Omega[-3])/2
            elif prof == 'rho':
                rho = prof_dict['rho']
                rho0 = reader.get_1d_list(scanparam = self.initparams[f'{reader}'], variable = 'rho0', spar_list = self.spar_lists[f'{reader}'][idx], paramSpecs = scan, _returnBoth = False)
                rho = rho*rho0

                step = (rho[0] - rho[-3])/2
                avg = (rho[0] + rho[-3])/2
            elif prof == 'beta':
                P = prof_dict['p']
                eps_a = 1/R0
                P = 2*P/eps_a**2

                step = (P[0] - P[-3])/2
                avg = (P[0] + P[-3])/2

            steps.append(step)
            avgs.append(avg)
        return steps, avgs

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

        # load x data
        if self['x_axis_type'] == 'peakedness':
            x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], 'peakedness', spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
        elif self['x_axis_type'] == 'peakedness_anal':
            x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], 'peakedness_anal', spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
        elif self['x_axis_type'] == 'xstep_norm':
            if type(reader) == AEmhd.Reader.AEread:
                x_step = reader.get_1d_list(self.initparams[f'{reader}'], f"{self.xkeys[f'{reader}']}_step", spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                x_avg = reader.get_1d_list(self.initparams[f'{reader}'], f"{self.xkeys[f'{reader}']}_avg", spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False) # i think this returns an array even if x_avg is constant?
            elif type(reader) == lunamhd.lunaReader.lunaRead:
                x_step = self.load_stepavg(reader, f"{self.xkeys[f'{reader}']}", scan = scan)[0]
                x_avg = self.load_stepavg(reader, f"{self.xkeys[f'{reader}']}", scan = scan)[1]
            x_vals = [i/j for i, j in zip(x_step, x_avg)]
        else:
            if self.xkey0s[f'{reader}'] == self.scankeys[f'{reader}']:
                x_vals = self.spar_lists[f'{reader}']
            else:
                x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.xkeys[f'{reader}'], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)

        # load y data
        if type(reader) == lunamhd.lunaReader.lunaRead:
            y_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
            # converting from gamma/wA0 to gamma/wAavg as used in analytic work
            if self.wA_avgNorm:
                rhosteps, rhoavgs = self.load_stepavg(reader, 'rho', scan = scan) # there is a rhostep and rhoavg value associated with every gam value
                conversion = [sqrt(j/(i+j)) for i,j in zip(rhosteps, rhoavgs)]
                y_vals = [i*j for i,j in zip(y_vals, conversion)]

            gam_vals = [i.real*10 for i in y_vals]
            a_gam_vals = None
            
            if self.scankeys[f'{reader}'] in ['beta0', 'beta1', 'beta_avg', 'beta_step'] and self['x_axis_type'] != 'xstep_norm':
                x_vals = [i*100 for i in x_vals]
            elif self.scankeys[f'{reader}'] == 'mach' and self['rot_axis_type'] in ['omega', 'Omega']: # Change x_vals from mach to omegahat if needed, probably broken but i am not fixing this rn
                _, x_vals = reader.get_1d_list(scanparam=self.scankeys[f'{reader}'], variable=self.xkeys[f'{reader}'], paramSpecs=scan)
            
        elif type(reader) == AEmhd.Reader.AEread:         
            if reader.info['scantype'] == 'full':
                y_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'][0], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                asy_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'][1], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
                gam_vals = [i.imag for i in y_vals] # need to change for general eps_a
            elif reader.info['scantype'] == 'asy':
                asy_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.ykeys[f'{reader}'][1], spar_list = self.spar_lists[f'{reader}'],  paramSpecs = scan, _returnBoth = False)
            a_gam_vals = [i.imag for i in asy_vals] # need to change for general eps_a
            # Change x_vals from omega to mach0 or mach1, pretty sure this is broken atm because of new omega options
            if self.scankeys[f'{reader}'] == 'omega' and self['rot_axis_type'] in ['mach0', 'mach1']: 
                _, x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.xkeys[f'{reader}'], paramSpecs = scan)

        # for txtfname, txtf in self.txts:
        #     self._load_txt(txtfile = txtf)

        return x_vals, gam_vals, a_gam_vals

    def plot_vals(self, reader = None, readeridx = 0, txtfname = None, csvfname = None, scan = {}):
        if len(self['own_ls']) > 0:
            lstyle = self['own_ls'][readeridx]
        else:
            lstyle = '-'

        if len(self['markfreq']) > 0:
            markfreq = self['markfreq'][readeridx]
        else:
            markfreq=1

        # truncate plots if needed
        if self['plotrange']['pstart'] is None:
            pstart = 0
        else:
            pstart = self['plotrange']['pstart']
        if self['plotrange']['pstop'] is None:
            pstop = len(self._load_data(reader = reader, scan = scan)[0])
        else:
            pstop = self['plotrange']['pstop']

        if reader:
            if type(reader) == lunamhd.lunaReader.lunaRead:
                x_vals, gam_vals, _ = self._load_data(reader = reader, scan = scan)
                # if self.scanlabel == 'VENUS-MHD (2,1)':
                #     self.ax.plot(x_vals[:-1], gam_vals[:-1], lstyle, label=f'{self.scanlabel}', markersize=self['markersize'], markevery=markfreq)
                self.ax.plot(x_vals[pstart:pstop], gam_vals[pstart:pstop], lstyle, label=f'{self.scanlabel}', markersize=self['markersize'], markevery=markfreq)
            elif type(reader) == AEmhd.Reader.AEread:
                if self['AE_visible']['gam']:
                    x_vals, gam_vals, _ = self._load_data(reader = reader, scan = scan)
                    self.ax.plot(x_vals[pstart:pstop], gam_vals[pstart:pstop], lstyle, label=f'{self.scanlabel}', markersize=self['markersize'], markevery=markfreq)
                if self['AE_visible']['a_gam']:
                    x_vals, _, a_gam_vals = self._load_data(reader = reader, scan = scan)
                    self.ax.plot(x_vals[pstart:pstop], a_gam_vals[pstart:pstop], lstyle, label=f'{self.scanlabel}', markersize=self['markersize'])
        elif txtfname:
            x_vals, y_vals = loadtxt(self.txts[f'{txtfname}'], unpack = True)
            self.ax.plot(x_vals, y_vals, '--v', label = txtfname, markersize=self['markersize'])
        elif csvfname:
            x_vals, y_vals = loadtxt(self.csvs[f'{csvfname}'], delimiter = ',', unpack = True, skiprows=1)
            y_vals = [i*0.1 for i in y_vals] # need a more general way to renormalize this
            self.ax.plot(x_vals, y_vals, '--v', label = csvfname, markersize=self['markersize'])

        
    def draw_fig(self):
        for idx, reader in enumerate(self.readers):
            scans = self.scans[f'{reader}']
            if scans:
                for scan in scans:
                    self.plot_vals(reader = reader, readeridx = idx, scan = scan)
            else:
                self.plot_vals(reader = reader, readeridx = idx)

        for txtfname in self.txts.keys():
            self.plot_vals(txtfname = txtfname)

        for csvfname in self.csvs.keys():
            self.plot_vals(csvfname = csvfname)
                            
        self.ax.set_ylabel(self._y_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis'])
        self.ax.set_xlabel(self._x_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis']) 
        if self.scanlabel:
            self.ax.legend(loc=self['legend_loc'])
            self.ax.legend_.set_visible(self['visible']['legend'])
        
        self.fig.canvas.draw_idle()
