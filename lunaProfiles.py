import numpy as np

class lunaProf(object):
    def __init__self(mode_type = 'KH', profparams = {})
        self.mode_type = mode_type
        self.profparams = {'profile':'rho_t', 'rstep':0.5, 'drstep':0.15, 'n0':1, 'nu_n':2, 'qr':1., 'rs':0.3, 
            'q0':0.938, 'qs':None, 'nu_q':2.}

        def __getitem__(self, key):
        if key in self.profparams:
            return self.profparams[key]
        else:
            print(f"ERROR: {key} not found")

    def make_profile(x, y0 = 1, nu_y = 6, offset = 0.05, shape = 'step'):
        # shape can be step or curve?
        if shape == 'step':
            prof = .5*y0*(1 + np.tanh(self['rstep']**2 - x**2)/self['drstep']**2) + offset
        elif shape == 'decreasing_poly':
            prof = y0*(1 - x**nu_y)
        elif shape == 'increasing_poly':
            prof = 0 
        else:
            print('ERROR: invalid shape provided. Valid shapes are: step, decreasing_poly, increasing_poly')


# needs a thing that determines whether something is a spline or not
