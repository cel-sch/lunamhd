# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy

from matplotlib.pyplot import subplots, show, ion, axes, tight_layout
from matplotlib.widgets import Slider, Button

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenval', # ['eigenval','margin_stab']
                    'x_axis_type':'initparam',
                    'EV_visible':{'gam':True, 'wr':False, 'a_gam':True, 'a_wr':False},
                    'fig_type':'general', # ['paper', 'singleplot']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[8.5,6],
                                'paper':[4.5,3.5]},
                    'linestyles':{'plain':'-','asy':'D-'},
                    'markersize':2,
                    'visible':{'suptitle':False, 'legend':True}}

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
                
        self.labels = {'omega':'$\hat{Ω}$','beta':'$\hat{β}$','delq':'$\hat{Δq}$','rho':'ρ','sig1':'σ',
                       'eps_a':'$ε_a$','Gamma':'Γ','gam':'$\hat{γ}$','asygam':'$\hat{γ}_{asy}$',
                       'y_step':'$(y_0-y_1)/2$', 'y_avg':'$y_{avg}$','y0':'$y_0$', 'y1':'$y_1$',
                       'EV':'$\hat{ω}','a_EV':'$\hat{ω}_{asy}$','wr':'$\hat{w}_r$','asywr':'$\hat{w}_r _{asy}$'}
                
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
        
    def open_plot(self, plotfilename = None, ykey = None):
        # Creates figure and axes
        self.fig, self.ax = subplots(figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        self.fig.set_tight_layout(True)
             
        self.scans = self.reader.info['scans']
        if self['suptitle'] is None:
            suptitle = ''
            for key, val in self.reader.info['fixedparams'].items():
                suptitle += f"{self._getlabel(key)} = {val} "
            self.fig.suptitle(suptitle,fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        else:
            self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        
        self._load_x_axis(self['x_axis_type'])
        self._load_y_axis(self['y_axis_type'])
        
        ion()
        show()
        self.draw_fig()
        
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
            self.ykeys = ['EV', 'a_EV']
            self._y_ax_label = self._getlabel('gam')

    def plot_vals(self, scan = None):
        if scan: 
            self.scanlabel = [f'{self._getlabel(key)}={keyval}' for key, keyval in scan.items()] # empty for 1D scans
            self.scanlabel = ', '.join(self.scanlabel)
            self.asyscanlabel = None
        else:
            self.scanlabel = None
            scan = {} # to set paramSpecs to {}

        x_vals, y_vals = self.reader.get_1d_list(self.xkey, self.ykeys[0], paramSpecs = scan) # need to check what happens if paramSpecs = None
        _, asy_vals = self.reader.get_1d_list(self.xkey, self.ykeys[1], paramSpecs = scan)
            
        if self.reader.info['scantype'] == 'full':
            if self['EV_visible']['gam']:
                gam_vals = [i.imag for i in y_vals]
                self.ax.plot(x_vals, gam_vals, self.lstyle, label=f'DE {self.scanlabel}', markersize=self['markersize'])
                self.asy_lstyle = self['linestyles']['asy']
            if self['EV_visible']['wr']:
                wr_vals = [i.real for i in y_vals]
                self.ax.plot(x_vals, wr_vals, self.lstyle, label=f'DE {self.scanlabel}', markersize=self['markersize'])
            
        if self['EV_visible']['a_gam']:
            a_gam_vals = [i.imag for i in asy_vals]
            self.ax.plot(x_vals, a_gam_vals, self.asy_lstyle, label=f'AE {self.scanlabel}', markersize=self['markersize'])
        if self['EV_visible']['a_wr']:
            a_wr_vals = [i.real for i in asy_vals]
            self.ax.plot(x_vals, a_wr_vals, self.asy_lstyle, label=f'AE {self.scanlabel}', markersize=self['markersize'])
        
    def draw_fig(self):
        self.lstyle = self['linestyles']['plain']
        self.asy_lstyle = self.lstyle

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
