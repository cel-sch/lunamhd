# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap
from numpy import sqrt

from matplotlib.pyplot import subplots, show, ion, axes, tight_layout
from matplotlib.widgets import Slider, Button

import AEmhd
import lunamhd

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenval', # ['eigenval','margin_stab']
                    'x_axis_type':'initparam',
                    'rot_axis_type':'mach0', # ['mach0', 'mach1', 'omega', 'Omega', 'omegahat']
                    'AE_visible':{'gam':True, 'a_gam':False},
                    'fig_type':'general', # ['paper', 'general']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[8.5,6],
                                'paper':[4.5,3.5]},
                    'linestyles':{'plain':'-x','asy':'D-'},
                    'markersize':2,
                    'visible':{'suptitle':True, 'legend':True, 'grid':True}}

class plot_multi(object):
    def __init__(self, readers = [], settings = {}):
        self.readers = readers
        self.settings = {}
        defaults = deepcopy(default_settings)
        
        for key in defaults:
                self.settings[key] = defaults[key]        
        for key in settings:
            if key not in defaults:
                print(f'ERROR: {key} not found.')
            else:
                self.settings[key] = settings[key]
                
        self.labels = {'omega':'$\hat{Ω}$', 'omegahat':'$\hat{Ω}$', 'mach':'$\mathcal{M}_0$',
                        'mach0':'$\mathcal{M}_0$', 'mach1':'$\mathcal{M}_1$',
                        'beta0':'$β_0$','beta':'$β_0$','rmaj':'$R_0$','b0':'$B_0$', 'eps_a':'$ε_a$',
                        'rho':'ρ',
                        'asygam':'$\hat{γ}_{asy}$', 'asywr':'$\hat{w}_r _{asy}$', 'a_EV':'$\hat{ω}_{asy}$',
                        'gam':'$\hat{γ}$','wr':'$\hat{w}_r$','EV':'$\hat{ω}',
                        'drstep':'$Δr_{step}$', 'q0':'$q_0$', 'rationalm':'m'}
                        # beta is for loading from AEs, is unnormalized here (eps_a normalizations are all removed for comparison to VENUS)
                
        self.xkeys = {}
        self.ykeys = {}
        self.scankeys = {}

        self.open_plot()
                
    def __getitem__(self, key):
        if key in self.settings:
            return self.settings[key]
        else:
            print(f"ERROR: {key} not found")
            
    def _getlabel(self, varkey):
        def _get_drive_label(self, format, varkey):
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
                self._get_drive_label(format=varkey, varkey='ρ')
            elif profile in ['PT']:
                self._get_drive_label(format=varkey, varkey='T') # or beta?
            elif profile in ['rot']:
                self._get_drive_label(fomrat=varkey, varkey='$\hat{Ω}$')
        elif varkey in self.labels.keys():
            label = self.labels[varkey]
        else:
            label = varkey
        return label
 
    def save_plot(self, plotfilename = None):
        if plotfilename is None:
            plotfilename = f"{self.reader.filename}"
            plotfilename = ''.join([i for i in str(plotfilename)[:-4]])
        self.fig.savefig(plotfilename)
        
    def open_plot(self, ykey = None):
        # Creates figure and axes
        self.fig, self.ax = subplots(figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        self.fig.set_tight_layout(True)
             
        self.scans = {}
        for reader in self.readers:
            self.scans[f'{reader}'] = reader.info['scans']
            self._load_x_axis(reader, self['x_axis_type'], self['rot_axis_type'])
            self._load_y_axis(reader, self['y_axis_type'])
        self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])

        if self['visible']['grid']:
            self.ax.grid()
        
        ion()
        show()
        self.draw_fig()

    def _make_point_info(self, reader):
        info = ''
        for key, val in reader.info['fixedparams'].items():
            info += f"{self._getlabel(key)} = {val} "
        info = "\n".join(wrap(info, 50))
        return info
        
    def _load_x_axis(self, reader, axis_type, rot_axis_type):
        if axis_type not in ['initparam']:
            print("ERROR: axis_type not found, valid types ['initparam']")
            return
        if axis_type == 'initparam':
            scanparam = reader.info['scanorder'][0]
            if scanparam in ['mach', 'Omega', 'omega']:
                xkey = self._load_rot_axis(reader, rot_axis_type)
            else:
                xkey = scanparam
            self.xkeys[f'{reader}'] = xkey
            self.scankeys[f'{reader}'] = scanparam
            self._x_ax_label = self._getlabel(xkey) # for now: just take the xlabel from the last reader to read in, probably not a good system
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
            self._y_ax_label = self._getlabel('gam')

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

    def _load_data(self, reader, scan = {}):
        if scan:
            self.scanlabel = [reader.info['runid']]
            self.scanlabel.append([f'{self._getlabel(key)}={keyval}' for key, keyval in scan.items()]) # empty for 1D scans
            self.scanlabel = ', '.join(self.scanlabel) 
        else:
            self.scanlabel = [reader.info['runid']]
        self.lstyle = self['linestyles']['plain']
        self.asy_lstyle = self['linestyles']['asy']

        if type(reader) == lunamhd.lunaReader.lunaRead:
            x_vals, y_vals = reader.get_1d_list(scanparam=self.scankeys[f'{reader}'], variable=self.ykeys[f'{reader}'], paramSpecs=scan)
            gam_vals = [i.real for i in y_vals]
            a_gam_vals = None
            # Change x_vals from mach to omegahat if needed
            if self.scankeys[f'{reader}'] == 'mach' and self['rot_axis_type'] in ['omega', 'Omega']: 
                _, x_vals = reader.get_1d_list(scanparam=self.scankeys[f'{reader}'], variable=self.xkeys[f'{reader}'], paramSpecs=scan)
            
        elif type(reader) == AEmhd.Reader.AEread:
            if reader.info['scantype'] == 'full':
                x_vals, y_vals = reader.get_1d_list(scanparam=self.scankeys[f'{reader}'], variable=self.ykeys[f'{reader}'][0], paramSpecs=scan, _renormVENUS=True)
                gam_vals = [i.imag for i in y_vals]
                _, asy_vals = reader.get_1d_list(scanparam=self.scankeys[f'{reader}'], variable=self.ykeys[f'{reader}'][1], paramSpecs=scan, _renormVENUS=True)
                a_gam_vals = [i.imag for i in asy_vals]
            elif reader.info['scantype'] == 'asy':
                x_vals, asy_vals = reader.get_1d_list(scanparam=self.scankeys[f'{reader}'], variable=self.ykeys[f'{reader}'][1], paramSpecs=scan, _renormVENUS=True)
                a_gam_vals = [i.imag for i in asy_vals]
            # Change x_vals from omega to mach0 or mach1
            if self.scankeys[f'{reader}'] == 'omega' and self['rot_axis_type'] in ['mach0', 'mach1']: 
                _, x_vals = reader.get_1d_list(self.scankeys[f'{reader}'], self.xkeys[f'{reader}'], paramSpecs = scan)

        return x_vals, gam_vals, a_gam_vals

    def plot_vals(self, reader, scan = {}):
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
        
    def draw_fig(self):
        for reader in self.readers:
            scans = self.scans[f'{reader}']
            if scans:
                for scan in scans:
                    self.plot_vals(reader = reader, scan = scan)
            else:
                self.plot_vals(reader = reader)
                            
        self.ax.set_ylabel(self._y_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis'])
        self.ax.set_xlabel(self._x_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis']) 
        if self.scanlabel:
            self.ax.legend()
            self.ax.legend_.set_visible(self['visible']['legend'])
        
        self.fig.canvas.draw_idle()
