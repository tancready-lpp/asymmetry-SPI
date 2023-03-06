# This is the proper script where the analysis for step1 of my workplan will be made

import numpy as np
import time
from getdat import getdat8
from ppf import ppfdata
# In this way I can work on the storage in a file and then import the results
from pulses_datacollection import pulse, pseq_stef
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Here I define another nested dictionary which contains the radiation asymmetry info for each pulse.
# The dictionary is structured as follows:
# r_asy[pulse number] = {
# t_max: [float] = time instant at which POW Vertical is max
# RAF: [array] = Radiaton Asymmetry factor as function of powv_time
# rafmax: [float] = value of RAF at t_max
# quantity_tmax: [float] = value of 'quantity' at Ip_t_t_max
# quantity_t_tmax: [float] = closest instant of time to t_max for quantity_t timevector
#   }

beg = time.time()

r_asy = dict()
# print(pseq_stef)

for pnum in pseq_stef:
    r_asy[pnum] = dict()  # Radiation asymmetry dictionary for each pulse
    # Calculate the Radiation Asymmetry Factor (raf)
    powV = pulse[pnum]['powv']
    powH = pulse[pnum]['powh']
    raf = (powV - powH) / (powV + powH)
    err_raf = 0.2 * raf
    # find the time at which powv or powh is max for every pulse and use the time for the max value between the two pow
    Vmax = max(powV)
    Hmax = max(powH)
    tol = 1e-4
    if Vmax>Hmax:
        r_asy[pnum]['powv_tmax'] = Vmax
        for t_index, value in enumerate(powV):
            if abs(value - Vmax) < tol:
                t_max_ind = t_index
        tmax = pulse[pnum]['powv_t'][t_max_ind]
    else:
        r_asy[pnum]['powh_tmax'] = Hmax
        for t_index, value in enumerate(powH):
            if abs(value - Hmax) < tol:
                t_max_ind = t_index
        tmax = pulse[pnum]['powh_t'][t_max_ind]
    r_asy[pnum]['tmax'] = tmax  # time instant at which Bolometer pow (V or H) is max

    # calculate the raf at the t_max (index!) and store it under a new key of the dict
    raf_max = raf[t_max_ind]
    r_asy[pnum]['rafmax'] = raf_max
    r_asy[pnum]['err_rafmax'] = err_raf[t_max_ind]
    # now what I need to do is calculate the quantities in the pulse dict at the t_max for Bolo pow

    for key in pulse[pnum].keys():
        if "_t" not in key and key is not 'phi_n1':
            keymax = np.interp(tmax, pulse[pnum][key + '_t'], pulse[pnum][key])
            # keymax = interp1d()<<
            r_asy[pnum][key + '_tmax'] = keymax
            # r_asy[pnum][key + '_tmax'] = pulse[pnum][key+ '_t'][t_max_ind]


# try to calculate the mode phase by myself
# for n in [96239,96109]:
#     sinmax = np.interp(tmax, getdat8("C2-MHDF", n)[1], getdat8("C2-MHDF", n)[0] )
#     cosmax = np.interp(tmax, getdat8("C2-MHDG", n)[1], getdat8("C2-MHDG", n)[0] )
#     # r_asy[n]['phase_tmax'] = np.arctan(sinmax/cosmax)*(180./np.pi)
#     print(np.arctan(sinmax/cosmax)*(180./np.pi))#, pulse[n]['phi_n1'])
#     # print(pulse[n]['phase_tmax'] + pulse[n]['phi_n1'])





ran = np.random.choice(pseq_stef)
keys = r_asy[ran].keys()
raf = {key: [] for key in keys}  # used a dict comprehension assigning an empty list to every item in a given list

# In the following loop I store the _tmax values for every quantity and for every pulse in a dict. In this way a
# given key is a _tmax quantity name and its value is the list holding the actal _tmax values for every pulse so for
# every quantity I have all the pulses in a list; it is kind of the reverse dictionary that I used up to now
for key in keys:
    for p in pseq_stef:
        raf[key].append(r_asy[p][key])

raf['phi_n1'] = []
for p in pseq_stef:
    raf['phi_n1'].append(pulse[p]['phi_n1'])


plt.rcParams.update({'font.size':15})

plt.plot(raf['efcc15_tmax'], raf['phi_n1'], 'bo', label = 'efcc15')
plt.plot(raf['efcc37_tmax'], raf['phi_n1'], 'ro', label = 'efcc37')
plt.legend(loc= 0)
plt.title('OHMIC')
plt.xlabel('EFCCs current [kA]')
# plt.ylabel('Locked mode $\phi_{n=1}$ position respect to SPI [deg]')
plt.minorticks_on()

plt.show()



fin = time.time()
print('Computational Analysis time  is {:.3f}s'.format(fin - beg))

