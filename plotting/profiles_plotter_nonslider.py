# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap
from numpy import linspace, sqrt, array

from matplotlib.pyplot import subplots, show, ion, axes, tight_layout, rcParams, figure
from matplotlib.widgets import Slider, Button

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'profiles', 
                    'x_axis_type':'initparam',
                    'fig_type':'general', # ['paper', 'singleplot']
                    'fontsizes':{'general':{'title':14,'axis':12,'suptitle':20},
                                 'paper':{'title':10,'axis':9,'suptitle':12}},
                    'figsizes':{'general':[8.5,6],
                                'paper':[4.5,3.5]},
                    'linestyles':{'plain':'-x','asy':'D-'},
                    'markersize':2,
                    'visible':{'suptitle':True, 'title':True, 'legend':True, 'grid':True}}

class plot_profiles(object):
    def __init__(self, reader, varnrs, scanparam, spar_list, settings = {}):
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

        if scanparam is None:
            self.scanparam = self.reader.info['scanorder'][0]
        else: 
            self.scanparam = scanparam
        if spar_list is None:
            self.spar_list = self.reader.info['scanparams'][self.scanparam] # should only load in scan parameters which are not part of fixed parameter list
        else:
            self.spar_list = list(spar_list)
                
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
        # Creates figure and axes
        if self['visible']['grid']:
            rcParams['axes.grid'] = True

        self.figs = {}
        self.scans = self.reader.info['scans']
        if self.scans: # create multiple figs
            for scan in self.scans: 
                for i in self.spar_list:
                    EF_file, _ = self.reader.get_profiles_list(scanparam = self.xkey, spar_list = i, paramSpecs = scan, _returnBoth = False).items():
                    EF_file = EF_file.split('/')[-1]
                    self.figs[f'{EF_file}'] = figure(num=EF_file, figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
                    self.figs[f'{EF_File}'].set_tight_layout(True)

                    self._make_suptitle(self.figs[f'{EF_File}'])
        else:
            for i in self.spar_list:
                    EF_file = self.reader.get_profiles_list(scanparam = self.xkey, spar_list = i, paramSpecs = {}, _returnBoth = False).keys():
                    EF_file = EF_file.split('/')[-1]
                    self.figs[f'{EF_file}'] = figure(num=EF_file, figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
                    self.figs[f'{EF_File}'].set_tight_layout(True)

                    self._make_suptitle(self.figs[f'{EF_File}'])
        
        self._load_x_axis()
        self._load_y_axis(self['y_axis_type'])
        
        self.draw_fig()
        ion()
        show()

    def _make_suptitle(self, fig):
        if self['suptitle'] is None:
            suptitle = ''
                for key, val in self.reader.info['fixedparams'].items():
                    suptitle += f"{self._getlabel(key)} = {val} "
                if len(self.varnrs) == 1:
                    suptitle += f'Variable{self.varnrs[0]} '
                suptitle = "\n".join(wrap(suptitle, 50))
            settings['suptitle'] = suptitle

        fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        
    def _load_x_axis(self):
        self.xkey = self.scanparam # required for input to retrieve list of profile files
        
    def _load_y_axis(self, axis_type):
        if axis_type not in ['profile','profiles']:
            print("ERROR: axis_type not found, valid types: ['profiles','profile']")
            return
        if axis_type in ['profile','profiles']:
            self.ykeys = 'EF_file'

    def plot_vals(self, scan = None):
        for EF_file, prof_dict in self.reader.get_profiles_list(scanparam = self.xkey, spar_list = self.spar_list, paramSpecs = scan, _returnBoth = False).items():
            EF_file = EF_file.split('/')[-1]
            fig = self.figs['EF_file']

            s = prof_dict['s']
            R = prof_dict['R']
            R0 = prof_dict['R0']
            Z = prof_dict['Z']

            ax1 = fig.add_subplot(2,4,1)
            ax1.plot(s, prof_dict['temp'])
            ax1.set_ylabel(r'$T/T_0$')
            ax1.set_xlabel('s')

            ax2 = fig.add_subplot(2,4,2)
            ax2.plot(s, prof_dict['q'])
            ax2.set_ylabel(r'q')
            ax2.set_xlabel('s')

            ax3 = fig.add_subplot(2,4,3)
            ax3.plot(s, prof_dict['press'], label=r'$P(s)$')
            ax3.plot(s, prof_dict['press_rot'], label=r'$P(s)e^{U(R^2-R_0^2)}$')
            ax3.legend()
            ax3.tick_params('x', labelbottom=False)
            ax3.set_ylabel(r'$\bar{P}$ [Pa]')
            ax3.set_xlabel('s')

            ax4 = fig.add_subplot(2,4,4)
            ax4.plot(R[:,-1]*R0, Z[:,-1]*R0)
            ax4.plot(R[0,0]*R0, Z[0,0]*R0,'+')
            ax4.set_xlabel(r'R [m]')
            ax4.set_ylabel(r'Z [m]')
            ax4.set_aspect('equal')

            ax5 = fig.add_subplot(2,4,5)
            ax5.plot(s, prof_dict['U'])
            ax5.set_ylabel(r'$U$ $[m^{-2}]$')
            ax5.set_xlabel('s')

            ax6 = fig.add_subplot(2,4,6)
            ax6.plot(s, prof_dict['rho'])
            ax6.set_ylabel(r'$\bar{\rho}/\rho_0$')
            ax6.set_xlabel('s')

            ax7 = fig.add_subplot(2,4,7)
            ax7.plot(s, prof_dict['omega'])
            ax7.set_ylabel(r'$Ω/Ω_0$')
            ax7.set_xlabel('s')
            
    def draw_fig(self):
        if self.scans:
            for scan in self.scans: 
                print(scan)
                self.plot_vals(scan = scan)
        else:
            self.plot_vals()
        
        #self.fig.canvas.draw_idle()
