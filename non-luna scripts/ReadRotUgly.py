# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 10:08:30 2022

@author: celin
"""

import matplotlib.pyplot as plt
import numpy as np

path = "Results/IK/UglyScript"
filename = "Run4.txt"

datapathWin = f"{path}/{filename}"

################################ VENUS-MHD Data ###############################
eps_a = 0.1

vs, machs, gams = np.loadtxt(datapathWin, dtype=(np.cfloat)).transpose()

machs = [i.real for i in machs]
vs = [i.real for i in vs]
gams = [i.real for i in gams] # normalised to eps_a to match growth rate in paper

################################ Plotting #####################################

plt.figure()
plt.plot(vs[:], gams[:], label='VENUS-MHD')
#plt.xlabel("$\hat{Ω}$")
plt.xlabel("$v_0/v_a$")
plt.ylabel("$γ/ω_A$")
plt.title(f"{filename}")
plt.legend()
plt.show()

