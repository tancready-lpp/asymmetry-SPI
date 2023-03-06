import numpy as np
from ppf import ppfdata
import time
# from pulses_datacollection import pulse
from scipy.integrate import trapz
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


start = time.time()
plt.rcParams.update({'font.size':18})

shots = [96238,96231,96103,98088, 98092, 98093,98071]
user = 'md7913'
threshold = dict()
for p in shots:
    n = 10
    threshold[p] = dict()

    threshold[p]['raw'] = ppfdata(p,'MCRW', 'RAWS', uid=user)[0]
    threshold[p]['raw_t'] = ppfdata(p, 'MCRW', 'RAWS', uid=user)[2]
    threshold[p]['gaus'] = ppfdata(p, 'MCRW', 'FIT', uid=user)[0]

    peak_t = list(find_peaks(abs(threshold[p]['gaus']))[0]) #INDEX VALUE
    peaks = list(threshold[p]['raw'][peak_t])
    dt = []

    if len(peak_t)!=1:
        for l in range(len(peak_t)-1):
            dt.append((threshold[p]['raw_t'][peak_t[l+1]] - threshold[p]['raw_t'][peak_t[l]])*1e3)
            threshold[p]['delay'] = dt
    else:
        threshold[p]['delay'] = 'one peak'

    fw = ppfdata(p, 'MCRW', 'FWHM', uid=user)[0]
    # print(fw)
    if len(fw) != len(peak_t):
        fw = fw[:len(peak_t)]

    mass = []
    sigma = fw / 2.
    threshold[p]['sigma'] = sigma

    # --- Adiacent boundary method --- #
    k = []
    for l in range(len(peak_t)-1):
        if (n)*(sigma[l]+sigma[l+1])*1e3 > dt[l]:
            n=2*dt[l]/((sigma[l]+sigma[l+1])*1e3)
            k.append(n)

    for i, sig in enumerate(fw/2.):
        s_index = 0
        t_start = threshold[p]['raw_t'][peak_t[i]] - n*sig
        t_stop = threshold[p]['raw_t'][peak_t[i]] + n*sig
        while abs(threshold[p]['raw_t'][s_index] - t_start)>1e-5:
            s_index+=1
        e_index = s_index
        while abs(threshold[p]['raw_t'][e_index] - t_stop) > 1e-5:
            e_index += 1
        t_int = np.linspace(t_start, t_stop, e_index-s_index)
        integral = trapz(threshold[p]['gaus'][s_index:e_index], t_int,dx = 1e-8)
        mass.append(integral)

    # --- +-2 sigma method --- #
    # for i, peak in enumerate(peaks):
    #     sigma = fw / 2.
    #     threshold[p]['sigma'] = sigma
    #     n=2 #integrate in 2 std dev
    #     s_index =0
    #     t_start = threshold[p]['raw_t'][peak_t[i]] - n * threshold[p]['sigma'][i]
    #     t_stop = threshold[p]['raw_t'][peak_t[i]] + n*threshold[p]['sigma'][i]
    #     print(t_start)
    #     while abs(threshold[p]['raw_t'][s_index] - t_start) > 1e-5:
    #         s_index+=1
    #     e_index = s_index
    #     while abs(threshold[p]['raw_t'][e_index] - t_stop) > 1e-5:
    #         e_index += 1
    #     t_int = np.linspace(t_start, t_stop, e_index - s_index)
    #     integral = trapz(threshold[p]['raw'][s_index:e_index], t_int, dx=1e-8)
    #     mass.append(integral)

    totmass = sum(mass)
    rel_mass = mass/totmass
    threshold[p]['relmass'] = rel_mass
    print(p, '{:.2e} '.format(totmass), ' V*ms (gaus)','rel mass',rel_mass)



n=96103
plt.figure()
plt.suptitle(str(n))

plt.plot(threshold[n]['raw_t'] - 40., threshold[n]['raw'], 'b-', label = 'MCRW raw data')
plt.plot(threshold[n]['raw_t']-40., threshold[n]['gaus'],'r-', label = 'Gaussian fit from PPF')
# for i in range(len(peaks)):
#     plt.axvline(x=(threshold[n]['raw_t'][peak_t[i]] - 2*threshold[n]['sigma'][i]), c='orange')
#     plt.axvline(x=(threshold[n]['raw_t'][peak_t[i]] + 2*threshold[n]['sigma'][i]))
plt.xlabel('Time [s]')
plt.ylabel('V')

plt.legend(loc='lower right')
# plt.xlim(16.85,16.9)
plt.ylim(-2,1)
plt.show()

end= time.time()
print('Computational Threshold time is {:.3f}s'.format(end-start))