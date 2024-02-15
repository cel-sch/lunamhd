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

default_settings = {'suptitle': None,
                    'title':None,
                    'y_axis_type':'eigenval', # ['eigenval','margin_stab']
                    'x_axis_type':'initparam',
                    'rot_axis_type':'mach', # ['omega','mach']
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
    def __init__(self, readers = {}, settings = {}):
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
                
        self.labels = {'omega':'$\hat{Ω}$','mach':'$\mathcal{M}$','beta0':'$β_0$','rmaj':'$R_0$','b0':'$B_0$',
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
        
    def open_plot(self, plotfilename = None, ykey = None):
        # Creates figure and axes
        self.fig, self.ax = subplots(figsize=(self['figsizes'][f"{self['fig_type']}"][0],self['figsizes'][f"{self['fig_type']}"][1]))
        self.fig.set_tight_layout(True)
             
        self.scans = {}
        for readid, reader in self.readers.items:
            self.scans[readid] = reader.info['scans']
        self.fig.suptitle(self['suptitle'],fontsize=self['fontsizes'][f"{self['fig_type']}"]['suptitle'],visible=self['visible']['suptitle'])
        
        self._load_x_axis(self['x_axis_type'])
        self._load_y_axis(self['y_axis_type'])

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

    def _load_rot_axis(self, reader, scan = {}, rotvals, axis_type):
        if axis_type not in ['mach', 'omega', 'Omega']:
            print("ERROR: axis_type not found, valid types: ['mach', 'omega', 'Omega']")
            return
        if axis_type in ['omega', 'Omega']:
            if type(reader) == lunamhd.lunaReader.lunaRead:
                rotvals = rotvals # do things to rotvals to make them from M into omega
        elif axis_type in ['mach']:
            if type(reader) == AEmhd.Reader.AEread: # should take normalised omega
                #also need to double check validity of this depending on the profiles. maybe build this into the reader instead??? and what if beta is varying? should be covered by scan
                omega_sq = [i**2 for i in rotvals]
                mach_sq = omega_sq/reader('beta', paramSpecs=scan)
                rotvals = [sqrt(i) for i in mach_sq]
        return rotvals

    def _renorm_AE_data(self, reader, scan = {}, param, vals):
        # this may likely need editing as normalisations change
        # currently written for the case where EVs are normalised to some kind of wA in both AEs and VENUS
        if 'eps_a' in reader.info['scanorder']:
            # need to figure out what to do because then eps_a is varying
        else:
            if param in ['omega', 'Omega']:
                # change omega to mach (elsewhere?)
            elif param in ['delq', 'EV', 'a_EV', ]:
                vals = [i*reader('eps_a', paramSpecs=scan)]
            elif param in ['beta']:
                vals = [i*reader('eps_a', paramSpecs=scan)**2]

    def plot_vals(self, reader, scan = {}):
        if scan:
            self.scanlabel = [reader]
            self.scanlabel.append([f'{self._getlabel(key)}={keyval}' for key, keyval in scan.items()]) # empty for 1D scans
            self.scanlabel = ', '.join(self.scanlabel) 
        else:
            self.scanlabel = reader
        lstyle = self['linestyles']['plain']
        asy_lstyle = lstyle

        if type(reader) == lunamhd.lunaReader.lunaRead:
            x_vals, y_vals = self.reader.get_1d_list(self.xkey, self.ykeys[0], paramSpecs = scan)
            self.ax.plot(x_vals, y_vals, label='?', markersize=self['markersize'])
        elif type(reader) == AEmhd.Reader.AEread:
            x_vals, y_vals = self.readers[reader].get_1d_list(self.xkey, self.ykeys[0], paramSpecs = scan)
            _, asy_vals = self.readers[reader].get_1d_list(self.xkey, self.ykeys[1], paramSpecs = scan)
            if reader.info['scantype'] == ['full']:
                if self['AE_visible']['gam']:
                    gam_vals = [i.imag for i in y_vals]
                if self['AE_visible']['a_gam']:
                    a_gam_vals = [i.imag for i in asy_vals]
                    asy_lstyle = self['linestyles']['asy']
            if reader.info['scantype'] == ['asy']:
                a_gam_vals = [i.imag for i in asy_vals]
        #elif type(reader) == text file: read text file
        
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
