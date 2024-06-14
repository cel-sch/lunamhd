# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:33:53 2023

@author: celin
"""
from copy import deepcopy
from textwrap import wrap
from numpy import linspace, sqrt, array, asarray, float64

from matplotlib.pyplot import figure, axes, subplots, show, ion, tight_layout, rcParams
from matplotlib.widgets import Slider, Button
from matplotlib.ticker import FormatStrFormatter

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
                    'visible':{'suptitle':True, 'title':True, 'legend':True, 'grid':True, 'infoplot':True},
                    'sliders':True}

class plot_profiles(object):
    def __init__(self, reader, scanparam, spar_list, settings = {}):
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
            self.spar_list = self.reader.info['scanparams'][self.scanparam]
        else:
            self.spar_list = spar_list

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

        self.fig, self.ax = subplots(nrows=2, ncols=4, figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))

        self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        self.fig.set_tight_layout(True)

        self.fig2, self.ax2 = subplots(figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        self.fig2.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        self.fig2.set_tight_layout(True)

        self.scans = self.reader.info['scans']

        # Make slider
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

    def _make_point_info(self, scanval):
        info = ''
        for key, val in self.reader.info['fixedparams'].items():
            if type(val) in ['int', float, float64]:
                info += f"{self._getlabel(key)} = {val:.4f}\n"
            else:
                info += f"{self._getlabel(key)} = {val}\n"
        if type(scanval) in ['int', float, float64]:
            info += f"{self.scanparam} = {scanval:.4f}"
        else:
            info += f"{self.scanparam} = {scanval}"
        #info = "\n".join(wrap(info, 30))
        return info

    def _load_x_axis(self):
        self.xkey = self.scanparam # required for input to retrieve list of profile files
        
    def _load_y_axis(self, axis_type):
        if axis_type not in ['profile','profiles']:
            print("ERROR: axis_type not found, valid types: ['profiles','profile']")
            return
        if axis_type in ['profile','profiles']:
            self.ykeys = 'EF_file'

    def plot_vals(self, scan = {}):
        val = self.slider.val if self.slider is not None else 0
        for key, vals in self.reader.get_profiles_list(scanparam = self.xkey, spar_list = self.spar_list[val], paramSpecs = scan, _returnBoth = False).items():
            EF_file = key
            prof_dict = vals
        EF_file = EF_file.split('/')[-1]

        s = prof_dict['s']
        R = prof_dict['R']
        R0 = prof_dict['R0']
        Z = prof_dict['Z']

        self.ax[0,0].cla()
        self.ax[0,0].plot(s, prof_dict['temp'])
        self.ax[0,0].set_ylabel(r'$T/T_0$')
        self.ax[0,0].set_xlabel('s')

        self.ax[0,1].cla()
        self.ax[0,1].plot(s, prof_dict['q'])
        self.ax[0,1].set_ylabel(r'q')
        self.ax[0,1].set_xlabel('s')

        self.ax[0,2].cla()
        self.ax[0,2].plot(s, prof_dict['press'], label=r'$P(s)$')
        self.ax[0,2].plot(s, prof_dict['press_rot'], label=r'$P(s)e^{U(R^2-R_0^2)}$')
        self.ax[0,2].legend()
        self.ax[0,2].tick_params('x', labelbottom=False)
        self.ax[0,2].set_ylabel(r'$\bar{P}$ [Pa]')
        self.ax[0,2].set_xlabel('s')

        self.ax[0,3].cla()
        self.ax[0,3].plot(R[:,-1]*R0, Z[:,-1]*R0)
        self.ax[0,3].plot(R[0,0]*R0, Z[0,0]*R0,'+')
        self.ax[0,3].set_xlabel(r'R [m]')
        self.ax[0,3].set_ylabel(r'Z [m]')
        self.ax[0,3].set_aspect('equal')

        self.ax[1,0].cla()
        self.ax[1,0].plot(s, prof_dict['U'])
        self.ax[1,0].set_ylabel(r'$U$ $[m^{-2}]$')
        self.ax[1,0].set_xlabel('s')

        self.ax[1,1].cla()
        self.ax[1,1].plot(s, prof_dict['rho'])
        self.ax[1,1].set_ylabel(r'$\bar{\rho}/\rho_0$')
        self.ax[1,1].set_xlabel('s')

        self.ax[1,2].cla()
        self.ax[1,2].plot(s, prof_dict['omega'])
        self.ax[1,2].set_ylabel(r'$Ω/Ω_0$')
        self.ax[1,2].set_xlabel('s')

        self.ax[1,3].cla()
        self.ax[1,3].grid(False)
        self.ax[1,3].xaxis.set_visible(False)
        self.ax[1,3].yaxis.set_visible(False)
        self.ax[1,3].set_frame_on(False)
        self.ax[1,3].text(0,0,self._make_point_info(scanval = self.spar_list[val]))
        if not self['visible']['infoplot']:
            self.ax[1,3].set_axis_off()

        shaf_diff0 = self.reader.get_1d_list(scanparam = self.xkey, variable = 'shaf_diff0', spar_list = self.spar_list[val], paramSpecs = scan, _returnBoth = False)
        shaf_diff0 = shaf_diff0[0]

        self.ax2.cla()
        self.ax2.plot(s, shaf_diff0)
        self.ax2.set_ylabel(r'$\frac{|J^\phi\mathcal{J}-\Delta^*\psi_p(\mathcal{J}/R^2)|}{max[\Delta^*\psi_p(\mathcal{J}/R^2)]}$')
        self.ax2.set_xlabel('s')
            
    def draw_fig(self, slider_idx=None):
        if self.scans:
            for scan in self.scans: 
                self.plot_vals(scan = scan)
        else:
            self.plot_vals()
        
        self.fig.canvas.draw_idle()
        self.fig2.canvas.draw_idle()
