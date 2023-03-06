
#this script just neds to print at what time we have the spikes in BOLOMETER power
# so that I can coorelate them with fast camera signals from JUVIL


import numpy as np
import time
# In this way I can work on the storage in a file and then import the results
from pulses_datacollection import pulse, pseq_stef
from pulses_analysis import r_asy


tmaxV = []
tmaxH = []
for pnum in pseq_stef:
    # Calculate the Radiation Asymmetry Factor (raf)
    powV = pulse[pnum]['powv']
    powH = pulse[pnum]['powh']
    # find the time at which powv or powh is max for every pulse and use the time for the max value between the two pow
    Vmax = max(powV)
    Hmax = max(powH)
    tol = 1e-4
    #tmax_V
    r_asy[pnum]['powv_tmax'] = Vmax
    for t_index, value in enumerate(powV):
        if abs(value - Vmax) < tol:
            t_max_ind_V = t_index
    tmaxV.append(pulse[pnum]['powv_t'][t_max_ind_V])
    #tmax_H
    r_asy[pnum]['powh_tmax'] = Hmax
    for t_index, value in enumerate(powH):
        if abs(value - Hmax) < tol:
            t_max_ind_H = t_index
    tmaxH.append(pulse[pnum]['powv_t'][t_max_ind_H])

print(pseq_stef)
print('tmax VERT = {}'.format(tmaxV))
print('tmax HOR = {}'.format(tmaxH))
