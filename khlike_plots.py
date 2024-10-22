import numpy as np
import os
from AEmhd import *
from lunamhd import *

markersize = 4
cwd = os.getcwd()

# ### VENUS vs DE section, impurity plot
# read = lunaRead('imp_omstep_lowbeta2')
# read2 = AEread('imp_omstep_beta00.5')
# mp = read.multi_plot(readers = [read2], settings={'x_axis_type':'xstep_norm','rot_axis_type':'omega', 'fig_type':'paper', 'reader_labels':['VENUS-MHD', 'step model'], 'own_ls':['cD-','m-'], 'markersize':markersize})
# mp.save_plot('imp_omstep_VENUSvDE.png', cwd)

# ### VENUS vs DE section, all stepped plot
# read = lunaRead('as_omstep_lowbeta2')
# read2 = AEread('as_omstep_beta00.5')
# mp = read.multi_plot(readers = [read2], settings={'x_axis_type':'xstep_norm','rot_axis_type':'omega', 'fig_type':'paper', 'reader_labels':['VENUS-MHD', 'step model'], 'own_ls':['cD-','m-'], 'markersize':markersize, 'plotrange':{'pstart':1,'pstop':None}})
# mp.save_plot('as_omstep_VENUSvDE.png', cwd)

# ### param deps section, temp vs omega
# # omstep only vs omstep + tempstep
# read = lunaRead('tvo_rotonly_2_1_drstep0.17_highres')
# read2 = lunaRead('tvo_rhostep0_2_1_drstep0.17_highres')
# read3 = AEread('tvo_rotonly_2_1')
# read4 = AEread('tvo_rhostep0_2_1')
# mp = read.multi_plot(readers = [read2, read3, read4], settings={'x_axis_type':'xstep_norm','rot_axis_type':'omega', 'x_axis_lims':[0,0.55], 'y_axis_lims':[-0.5,9.5], 'fig_type':'paper', 'reader_labels':['VENUS $\hat{Ω}$ step only', 'VENUS $\hat{Ω}, T/T_0, \hat{β}$ step', 'analytic $\hat{Ω}$ step only', 'analytic $\hat{Ω}, T/T_0, \hat{β}$ step'], 'axis_labels':{}, 'own_ls':['cD-','mD-', 'c-', 'm-'], 'markersize':markersize, 'markfreq':[1,1,1,2]})
# mp.ax.plot(0.1,0,'bD')
# mp.save_plot('tvo_VENUS.png', cwd)

# # rhostep only
# read = lunaRead('tvo_rhostep_2_1_drstep0.17_highres')
# read2 = AEread('tvo_rhostep_2_1')
# mp = read.multi_plot(readers = [read2], settings={'x_axis_type':'xstep_norm','rot_axis_type':'omega', 'fig_type':'paper', 'reader_labels':['VENUS-MHD', 'step model'], 'own_ls':['bD-','b-'], 'markersize':markersize, 'markfreq':[2,1]})
# mp.save_plot('tvo_rhostep_VENUS.png', cwd)

# ### param deps section, density drive
# # omstep vs rhostep (1,1)
# read = lunaRead('as_om_rho0_avg1')
# read2 = lunaRead('as_rho_om0_avg1')
# read3 = AEread('as_om_rho0_avg1')
# read4 = AEread('as_rho_om0_avg1')
# mp = read.multi_plot(readers = [read2, read3, read4], settings={'x_axis_type':'xstep_norm', 'rot_axis_type':'omega', 'fig_type':'paper', 'reader_labels':['VENUS $\hat{γ}(\\tilde{Ω}_{step})$', 'VENUS $\hat{γ}(\\tilde{ρ}_{step})$', 'analytic $\hat{γ}(\\tilde{Ω}_{step})$', 'analytic $\hat{γ}(\\tilde{ρ}_{step})$'], 'axis_labels':{'x':'$\\tilde{Ω}_{step}, \\tilde{ρ}_{step}$'}, 'own_ls':['cD-','mD-', 'c-', 'm-'], 'markersize':markersize, 'markfreq':[1,1,1,1]})
# mp.save_plot('omstepvrhostep_q1.png', cwd)

# # omstep vs rhostep (2,1)
# read = lunaRead('as_om_rho0_avg1_2_1')
# read2 = lunaRead('as_rho_om0_avg1_2_1')
# read3 = AEread('as_om_rho0_avg1_2_1')
# read4 = AEread('as_rho_om0_avg1_2_1')
# mp = read.multi_plot(readers = [read2, read3, read4], settings={'x_axis_type':'xstep_norm', 'rot_axis_type':'omega', 'fig_type':'paper', 'reader_labels':['VENUS $\hat{γ}(\\tilde{Ω}_{step})$', 'VENUS $\hat{γ}(\\tilde{ρ}_{step})$', 'analytic $\hat{γ}(\\tilde{Ω}_{step})$', 'analytic $\hat{γ}(\\tilde{ρ}_{step})$'], 'axis_labels':{'x':'$\\tilde{Ω}_{step}, \\tilde{ρ}_{step}$'}, 'own_ls':['cD-','mD-', 'c-', 'm-'], 'markersize':markersize, 'markfreq':[1,1,1,1]})
# mp.save_plot('omstepvrhostep_q2.png', cwd)

# ### mode nr dependence
# # (n,n)
# read = lunaRead('as_omstep_1_1')
# read2 = lunaRead('as_omstep_2_2')
# read3 = lunaRead('as_omstep_3_3')
# read4 = AEread('as_omstep_1_1')
# read5 = AEread('as_omstep_2_2')
# read6 = AEread('as_omstep_3_3')
# mp = read.multi_plot(readers = [read2, read3, read4, read5, read6], settings={'x_axis_type':'xstep_norm', 'rot_axis_type':'omega', 'fig_type':'paper', 'reader_labels':['VENUS (1,1)', 'VENUS (2,2)', 'VENUS (3,3)', 'analytic (1,1)', 'analytic (2,2)', 'analytic (3,3)'], 'own_ls':['cD-','mD-', 'yD-', 'c-', 'm-', 'y-'], 'markersize':markersize, 'markfreq':[1,1,1,1,1,1]})
# mp.save_plot('as_omstep_q1.png', cwd)

# (2n, n)
# read = lunaRead('as_omstep_2_1')
# read2 = lunaRead('as_omstep_4_2')
# read3 = AEread('as_omstep_2_1')
# read4 = AEread('as_omstep_4_2')
# read5 = AEread('as_omstep_6_3')
# mp = read.multi_plot(readers = [read2, read3, read4, read5], settings={'x_axis_type':'xstep_norm', 'rot_axis_type':'omega', 'fig_type':'paper', 'x_axis_lims':[-0.5,0.85], 'y_axis_lims':[-0.5,1], 'reader_labels':['VENUS (2,1)', 'VENUS (4,2)', 'analytic (2,1)', 'analytic (4,2)', 'analytic (6,3)'], 'own_ls':['cD-','mD-', 'c-', 'm-', 'y-'], 'markersize':markersize, 'markfreq':[1,1,1,1,1]})
# mp.save_plot('as_omstep_q2.png', cwd)

### beta dependence
# beta_avg scan
read = lunaRead('as_betaavg2')
read2 = AEread('as_betaavg')
mp = read.multi_plot(readers = [read2], settings={'x_axis_type':'initparam','rot_axis_type':'omega', 'fig_type':'paper', 'axis_labels':{'x':'$\hat{β}_{avg}$'}, 'reader_labels':['VENUS-MHD', 'step model'], 'own_ls':['cD-','m-'], 'markersize':markersize, 'markfreq':[1,1]})
mp.save_plot('as_betaavg.png', cwd)

# # beta_step scan
read = lunaRead('as_betastep_1_12')
read2 = lunaRead('as_betastep_2_1_long')
read3 = AEread('as_betastep_1_12')
read4 = AEread('as_betastep_2_1')
mp = read.multi_plot(readers = [read2, read3, read4], settings={'x_axis_type':'initparam','rot_axis_type':'omega', 'fig_type':'paper', 'reader_labels':['VENUS-MHD (1,1)', 'VENUS-MHD (2,1)', 'step model (1,1)', 'step model (2,1)'], 'own_ls':['cD-','mD-', 'c-', 'm-'], 'markersize':markersize, 'markfreq':[1,1,1,1], 'plotrange':{'pstart':0,'pstop':-1},'legend_loc':'upper left'})
mp.save_plot('as_betastep.png', cwd)