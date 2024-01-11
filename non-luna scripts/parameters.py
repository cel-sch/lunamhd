# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 13:22:39 2023

@author: celin
"""

import matplotlib.pyplot as plt
import numpy as np

### STEPPED ROTATION ###
s = np.linspace(0,1,100)
s2 = s**2

rstep = 0.5
drstep = 0.4

# rotation profile
Omega = 1.5*(1 + np.tanh((rstep**2 - s2)/drstep**2))

# density and pressure profiles are based on what i've actually implemented in venus-mhd, not on the 2013 ppcf ones
#n0 = 2
#n_ = n0*(1.-s**6)
#p = n_

# q profile
q0 = 1.01
qr = 1
qs = 2
nu_q = 8

q = q0+qs*s**nu_q

### STEPPED RHO_P ###

# rotation profile
Omega = 2*(1.-s**6)

# density and pressure profiles
n_ = 1.5*(1 + np.tanh((rstep**2 - s2)/drstep**2))
p = n_

# plot
fig, ax = plt.subplots()
ax.plot(s,Omega,label='$Ω(r)/Ω(0)$')
#ax.plot(s,n_,'-D',color='orange',label='$ρ(r)/ρ(0)$',markevery=6)
ax.plot(s,p,color='firebrick',label='$p(r)/p(0)$')
ax.plot(s,q,color='orange',label='$q(r)$')
ax.set_xlabel('r',fontsize=12)
ax.set_ylabel('',fontsize=12)
plt.text(0.08,2.6,'$ρ(r)/ρ(0)$',fontsize=12)
#plt.text(0.08,2.2,'$Δr_0$=0.4',fontsize=12)
plt.text(0.18,1.8,'$Ω(r)/Ω(0)$',fontsize=12)
plt.text(0.08,2.4,'$p(r)/p(0)$',fontsize=12)
plt.text(0.9,1.6,'$q(r)$',fontsize=12)
plt.text(0.52,2.7,'$r_0$',fontsize=12)
ax.axvline(0.5, ymin=0, ymax=1, linestyle='--',color='black')
fig.set_size_inches(7, 5.5)
#plt.legend()