# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap

from matplotlib.pyplot import subplots, show, ion, axes, tight_layout
from matplotlib.widgets import Slider, Button

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenval', # ['eigenval','margin_stab']
                    'x_axis_type':'initparam',
                    'EV_visible':{'gam':True, 'wr':False},
                    'EV_guess':True,
                    'fig_type':'general', # ['paper', 'singleplot']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[8.5,6],
                                'paper':[4.5,3.5]},
                    'linestyles':{'plain':'-x','asy':'D-'},
                    'markersize':2,
                    'visible':{'suptitle':True, 'legend':True, 'grid':True}}

class plot_growth(object):
    def __init__(self, reader, settings = {}):
        self.reader = reader
        self.settings = {}
        defaults = deepcopy(default_settings)
        
        for key in defaults:
                self.settings[key] = defaults[key]        
        for key in settings:
            if key not in defaults:
                print(f'ERROR: {key} not found.')
            else:
                self.settings[key] = settings[key]
                
        self.labels = {'mach':'$\mathcal{M}$','beta0':'$β_0$','rmaj':'$R_0$','b0':'$B_0$',
                       'gam':'$\hat{γ}$','wr':'$\hat{w}_r$','EV':'$\hat{ω}',
                       'drstep':'$Δr_{step}$', 'q0':'$q_0$', 'rationalm':'m'}
                
        self.xkey = None
        self.ykeys = None

        self.open_plot()
                
    def __getitem__(self, key):
        if key in self.settings:
            return self.settings[key]
        else:
            print(f"ERROR: {key} not found")
            
    def _getlabel(self, var):
        try:
            label = self.labels[var]
        except:
            label = var
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
             
        self.scans = self.reader.info['scans']
        if self['suptitle'] is None:
            suptitle = self._make_point_info()
            self.fig.suptitle(suptitle,fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        else:
            self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        
        self._load_x_axis(self['x_axis_type'])
        self._load_y_axis(self['y_axis_type'])

        if self['visible']['grid']:
            self.ax.grid()
        
        ion()
        show()
        self.draw_fig()

    def _make_point_info(self):
        info = ''
        for key, val in self.reader.info['fixedparams'].items():
            info += f"{self._getlabel(key)} = {val} "
        info = "\n".join(wrap(info, 50))
        return info
        
    def _load_x_axis(self, axis_type):
        if axis_type not in ['initparam']:
            print("ERROR: axis_type not found, valid types ['initparam']")
            return
        if axis_type == 'initparam':
            self.xkey = self.reader.info['scanorder'][0]
        self._x_ax_label = self._getlabel(self.xkey)
        
    def _load_y_axis(self, axis_type):
        if axis_type not in ['eigenval']:
            print("ERROR: axis_type not found, valid types: ['eigenval']")
            return
        if axis_type == 'eigenval':
            self.ykeys = 'EV'
            self._y_ax_label = self._getlabel('gam')

    def plot_vals(self, scan = None):
        if scan: 
            self.scanlabel = [f'{self._getlabel(key)}={keyval}' for key, keyval in scan.items()] # empty for 1D scans
            self.scanlabel = ', '.join(self.scanlabel)
        else:
            self.scanlabel = None
            scan = {} # to set paramSpecs to {}

        x_vals, y_vals = self.reader.get_1d_list(self.xkey, self.ykeys, paramSpecs = scan) # need to check what happens if paramSpecs = None
            
        if self['EV_visible']['gam']:
            gam_vals = [i.real for i in y_vals]
            self.ax.plot(x_vals, gam_vals, self.lstyle, label=f'{self.scanlabel}', markersize=self['markersize'])
            if 'EV_guess':
                _, gam_guess_vals = self.reader.get_1d_list(self.xkey, 'EVguess', paramSpecs = scan)
                gam_guess_vals = [i.real for i in gam_guess_vals]
                self.ax.plot(x_vals, gam_guess_vals, self.lstyle, label=f'{self.scanlabel} guess', markersize=self['markersize'])
        if self['EV_visible']['wr']:
            wr_vals = [i.imag for i in y_vals]
            self.ax.plot(x_vals, wr_vals, self.lstyle, label=f'{self.scanlabel}', markersize=self['markersize'])
        
    def draw_fig(self):
        self.lstyle = self['linestyles']['plain']

        if self.scans:
            for scan in self.scans:
                print(scan)
                self.plot_vals(scan = scan)
        else:
            self.plot_vals()
                    
        self.ax.set_ylabel(self._y_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis'])
        self.ax.set_xlabel(self._x_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis']) 
        if self.scanlabel:
            self.ax.legend()
            self.ax.legend_.set_visible(self['visible']['legend'])
        
        self.fig.canvas.draw_idle()
