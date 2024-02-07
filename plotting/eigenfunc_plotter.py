# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap
from numpy import linspace

from matplotlib.pyplot import subplots, show, ion, axes, tight_layout
from matplotlib.widgets import Slider, Button

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenfunc', # ['eigenval','margin_stab']
                    'x_axis_type':'initparam',
                    'fig_type':'general', # ['paper', 'singleplot']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[8.5,6],
                                'paper':[4.5,3.5]},
                    'linestyles':{'plain':'-x','asy':'D-'},
                    'markersize':2,
                    'visible':{'suptitle':True, 'title':True, 'legend':True, 'grid':True}}

class plot_EF(object):
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
        
    def open_plot(self, subplot_nr, plotfilename = None):
        # Creates figure and axes
        self.subplot_nr = subplot_nr

        def row_col(subplot_nr):
            # Finds the number of rows and columns best to use for EF subplots depending on the number of subplots
            def close_factors(subplot_nr):
                # Find closest pair of factors for subplot_nr to have a nice arrangement of the subplots
                factor1 = 0
                factor2 = subplot_nr
                while factor1+1 <= factor2:
                    factor1 += 1
                    if number % factor1 == 0:
                        factor2 = number // factor1
                return factor1, factor2

            while True:
                factor1, factor2 = close_factors(subplot_nr)
                if 0.5*factor1 <= factor2:
                    break
                subplot_nr += 1
            return factor1, factor2

        rownr, colnr = row_col(subplot_nr = subplot_nr)
        self.fig, self.axs = subplots(nrows = rownr, ncols = colnr, figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        self.fig.set_tight_layout(True)
             
        self.scans = self.reader.info['scans']
        if self['suptitle'] is None:
            suptitle = ''
            for key, val in self.reader.info['fixedparams'].items():
                suptitle += f"{self._getlabel(key)} = {val} "
            suptitle = "\n".join(wrap(suptitle, 50))
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
        
    def _load_x_axis(self, axis_type):
        if axis_type == 'initparam':
            self.xkey = self.reader.info['scanorder'][0] # required for input to retrieve list of eigenfuncs
        self.r = linspace(0.,1.,10000)
        self._x_ax_label = 'r'
        
    def _load_y_axis(self, axis_type):
        if axis_type not in ['eigenfunc']:
            print("ERROR: axis_type not found, valid types: ['eigenfunc']")
            return
        if axis_type == 'eigenfunc':
            self.ykeys = 'EF_file'
            self._y_ax_label = 'a.u.'

    def plot_vals(self, varnrs, scan = None):
        if scan: 
            self.scanlabel = [f'{self._getlabel(key)}={keyval}' for key, keyval in scan.items()] # empty for 1D scans
            self.scanlabel = ', '.join(self.scanlabel)
        else:
            self.scanlabel = None
            scan = {} # to set paramSpecs to {}

        EF_keys, EF_ms = self.reader.get_eigenfunc_list(varnrs = varnrs, self.xkey, self.ykeys, paramSpecs = scan) # need to check what happens if paramSpecs = None
        
        for i, EF_key in enumerate(EF_keys): # subplot_nr should be the same as the length of EF_keys, print error if not the case
            for m_val, EF in EF_ms.items():
                mode = EF
                vars()['ax'+f'{i}'].plot(self.r, mode, label=f'{m_val}') 
                if self['visible']['title']:
                    vars()['ax'+f'{i}'].set_title(f'{EF_key}') 
        
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
