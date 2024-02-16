# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap
from numpy import linspace, sqrt, array

from matplotlib.pyplot import figure, axes, subplots, show, ion, tight_layout, rcParams
from matplotlib.widgets import Slider, Button

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenfunc', 
                    'x_axis_type':'initparam',
                    'fig_type':'general', # ['paper', 'singleplot']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[8.5,6],
                                'paper':[4.5,3.5]},
                    'linestyles':{'plain':'-x','asy':'D-'},
                    'markersize':2,
                    'visible':{'suptitle':True, 'title':True, 'legend':True, 'grid':True},
                    'sliders':True}

class plot_EF(object):
    def __init__(self, reader, varnrs, scanparam, spar_list, settings = {}):
        self.reader = reader
        self.settings = {}
        self.slider = None
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

        if scanparam is None:
            self.scanparam = self.reader.info['scanorder'][0]
        else: 
            self.scanparam = scanparam
        if spar_list is None:
            self.spar_list = self.reader.info['scanparams'][self.scanparam] # should only load in scan parameters which are not part of fixed parameter list
        else:
            self.spar_list = spar_list
        self.varnrs = self.reader._make_list(varnrs)

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
        
    def open_plot(self):
        # Creates figure, axes and slider
        if self['visible']['grid']:
            rcParams['axes.grid'] = True

        def row_col(subplot_nr):
            # Finds the number of rows and columns best to use for EF subplots depending on the number of subplots
            if subplot_nr == 1:
                return 1, 1

            def close_factors(subplot_nr):
                # Find closest pair of factors for subplot_nr to have a nice arrangement of the subplots
                factor1 = 0
                factor2 = subplot_nr
                while factor1+1 <= factor2:
                    factor1 += 1
                    if subplot_nr % factor1 == 0:
                        factor2 = subplot_nr // factor1
                return factor1, factor2

            while True:
                factor1, factor2 = close_factors(subplot_nr)
                if 0.5*factor1 <= factor2:
                    break
                subplot_nr += 1
            return factor1, factor2

        self.rownr, self.colnr = row_col(subplot_nr = len(self.varnrs))
        self.fig, self.ax = subplots(nrows=self.rownr, ncols=self.colnr, figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        if self['suptitle'] is None:
            self.settings['suptitle'] = self._make_point_info()
        self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        self.fig.set_tight_layout(True)
        self.ax = array(self.ax)
        self.ax = self.ax.flatten()
        
        # Make slider    
        self.scans = self.reader.info['scans']
        if self['sliders']:
            self.slider_fig, self.slider_ax = subplots()
            self.slider_ax.grid(False)
            self.slider_ax.xaxis.set_visible(False)
            self.slider_ax.yaxis.set_visible(False)
            self.slider_ax.set_frame_on(False)
            slider_subax = self.slider_ax.axes.add_child_axes(axes([0.1, 0.47, 0.8, 0.06]))
            if self.scans:
                self.slider = Slider(slider_subax, label='point', valmin=0, valmax=len(self.scans)*len(self.spar_list)-1, valstep=1)
                self.slider.on_changed(self.draw_fig)
            else:
                self.slider = Slider(slider_subax, label='point', valmin=0, valmax=len(self.spar_list)-1, valstep=1)
                self.slider.on_changed(self.draw_fig)
        
        self._load_x_axis()
        self._load_y_axis(self['y_axis_type'])
        
        ion()
        show()
        self.draw_fig()

    def _make_point_info(self):
        # Does not currently include value of scanparam
        info = ''
        for key, val in self.reader.info['fixedparams'].items():
            info += f"{self._getlabel(key)} = {val}\n"
        #info += f"{self.scanparam} = {scanval}"
        info = "\n".join(wrap(info, 50))
        return info
        
    def _load_x_axis(self):
        self.xkey = self.scanparam # required for input to retrieve list of eigenfuncs
        self.r = linspace(0.,1.,102) # need to read this properly from the h5 file probably
        self.r = sqrt(self.r)
        self._x_ax_label = 'r'
        
    def _load_y_axis(self, axis_type):
        if axis_type not in ['eigenfunc']:
            print("ERROR: axis_type not found, valid types: ['eigenfunc']")
            return
        if axis_type == 'eigenfunc':
            self.ykeys = 'EF_file'
            self._y_ax_label = 'a.u.'

    def plot_vals(self, scan = {}):
        val = self.slider.val if self.slider is not None else 0
        for key, vals in self.reader.get_eigenfunc_list(varnrs = self.varnrs, scanparam = self.xkey, spar_list = self.spar_list[val], paramSpecs = scan, _returnBoth = False).items():
            EF_file = key
            EF = vals # dictionary like {'var_1':..., 'var_2':..., ...}

        for i, (var_key, var_allms) in enumerate(EF.items()):
            self.ax[i].cla()
            for m_val, mode in var_allms.items():
                self.ax[i].plot(self.r, mode, label=f'{m_val}') 
                if self['visible']['title']:
                    self.ax[i].set_title(f'Variable {self.varnrs[i]+1}, {self._getlabel(self.xkey)} = {self.spar_list[val]}')
                self.ax[i].legend()
                self.ax[i].legend_.set_visible(self['visible']['legend'])
                
                self.ax[i].set_ylabel(self._y_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis'])
                self.ax[i].set_xlabel(self._x_ax_label,fontsize=self['fontsizes'][f"{self['fig_type']}"]['axis'])
            
    def draw_fig(self, slider_idx=None):
        if self.scans:
            for scan in self.scans:
                print(scan)
                self.plot_vals(scan = scan)
        else:
            self.plot_vals()
        
        self.fig.canvas.draw_idle()
