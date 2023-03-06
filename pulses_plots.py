
from pulses_datacollection import  pseq_stef  # used for pulse info plot
from pulses_analysis import raf
from pulses_pellets import pellet
from q2 import pulse
import numpy as np  # used for pulse info plot
import matplotlib.pyplot as plt
import time
from ppf import ppfdata
from scipy.signal import find_peaks
# from broken_threshold import threshold
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



plt.rcParams.update({'font.size':15})
s = time.time()
# n = np.random.choice(pseq_stef)
n=98071
user = 'md7913'
#raf.keys() = [Ip_tmax, nei3_tmax, nef3_tmax, efcc15_tmax, efcc37_tmax, loca_tmax, phi_n1, phase_tmax, powh_tmax, powv_tmax]

xlab = ['Plasma Current [MA] at tmax',  # 0
        'Toroidal Field [T] at tmax', # 1
        'Integrated density from LID3 [10e20 m^-2] at tmax',  # 2
        'Integrated Density from Polarimeter from LID3 [10e20 m^-2] at tmax',  # 3
        'Efcc current [kA] at tmax',  # 4
        'Locked mode n=1 Amplitude [mT] at tmax',  # 5
        'Phi_n=1 locked mode O-point phase [deg] respect to SPI angular position',  # 6
        'Pow Vertical bolometer [GW] at tmax',  # 7
        'Pow Horizontal bolometer [GW] at tmax'  # 8
        ]


# for p in pseq_stef:
#         print('pulse ' + str(p) + ': phi = ' + str(pulse[p]['phi_n1']))
#
# plt.figure(figsize=(16,9), facecolor='white')
# # plot of ALL pulses for RAFMAX
# # plt.plot(raf['efcc15_tmax'], raf['rafmax'], 'bo', label = 'Efcc15 current')
# # # plt.plot(raf['efcc37_tmax'], raf['rafmax'], 'ro', label = 'Efcc37 current')
# # plt.xticks(np.arange(-90,360,45))
#
# plt.plot(raf['phi_n1'][0:2], raf['rafmax'][0:2], 'bo', label = '100%Ne w/out EFCC', markersize= 12)
# plt.plot(raf['phi_n1'][0:12], raf['rafmax'][0:12], 'bD', label = '100%Ne',markersize= 12)
# plt.plot(raf['phi_n1'][12:15], raf['rafmax'][12:15], 'gD', label = '100%Ne broken',markersize= 12)
# plt.plot(raf['phi_n1'][15:], raf['rafmax'][15:], 'rs', label = '10%Ne',markersize= 12)
#
# # plt.plot(raf['phase_tmax'][0:2], raf['rafmax'][0:2], 'bo', label = '100%Ne w/out EFCC', markersize= 12)
# # plt.plot(raf['phase_tmax'][2:12], raf['rafmax'][2:12], 'bD', label = '100%Ne',markersize= 12)
# # plt.plot(raf['phase_tmax'][12:15], raf['rafmax'][12:15], 'gD', label = '100%Ne broken',markersize= 12)
# # plt.plot(raf['phase_tmax'][15:], raf['rafmax'][15:], 'rs', label = '10%Ne',markersize= 12)

# # plt.ylabel('Radiation Asymmetry Fraction at tmax')
# plt.xlabel(xlab[6])
# plt.legend(loc = 'best', numpoints = 1 , fontsize = 15, markerscale = 0.75)
#
# plt.xlim(-100, 250)

# to annonate the pulse number next to the data point
# for i, p in enumerate(pseq_stef):
#         plt.annotate(str(int(p)), (raf['phi_n1'][i], raf['rafmax'][i]), textcoords = 'offset points', xytext =(5,5))

q = ['Ip','Te92','nef3','efcc','loca']#'nef3',',]#'my_phi_n1']#,'Bt',c 'powv','powh'
# q = []# 'powv', 'powh']' #oss: for density you can use nef3 or nei3 'Te_q2' 'Bt' ,
u = ['[MA]','[keV]','$[10^{20}m^{-2}$]','[kA]','[mT]']#,]#'[deg]']#'[$10^{20}m^{-2}$]','[deg], '[T]','[keV]','[GW]', '[GW]'
# u =[]



mc = abs(ppfdata(n, 'MCRW', 'FIT', uid=user)[0])
mc_t = ppfdata(n, 'MCRW', 'RAWS', uid=user)[2]
peak_t = find_peaks(mc)[0] #INDEX VALUE
# print(peak_t)
mc_t_pass = mc_t[peak_t]
# print(mc_t_pass)

minor = [.5,.2,.1,.1,.5]

efcc_numb = ['15']#,'37']
efcc_labels = ['efcc15']#  , 'efcc37']

# # plot of a random pulse of common quantities
fig, ax = plt.subplots(len(q), sharex=True, figsize=(12,12), facecolor='white')
fig.suptitle('PULSE ' + str(n) )#+'\nq = 2 ch is '+ str(pulse[n]['q2_kk3_ch']))
lab = ['Ip', 'Te core', 'ne core','powV','powH']# 'efcc','Br n=1'


for i, quant in enumerate(q):
    if quant is not 'efcc':
        ax[i].plot(pulse[n][quant + '_t']-40., pulse[n][quant], linewidth = 1)
        ax[i].set_ylabel(lab[i] + '\n' + u[i])
        ax[i].xaxis.set_minor_locator(MultipleLocator(.1))
        ax[i].yaxis.set_minor_locator(MultipleLocator(minor[i]))
        ax[i].axvline(x=mc_t_pass-40., c='r')
        ax[i].axvline((pellet[n]['plasma_t'] - 40.), c='r', ls='--', lw=1)
        # start, end = ax[i].get_xlim()
        # ax[i].set_xticks(np.linspace(start, end, 1))
    else:
        for j, number in enumerate(efcc_numb):
            ax[i].plot(pulse[n][quant + number + '_t']-40., pulse[n][quant+number], label = efcc_labels[j], linewidth = 2)
            ax[i].set_ylabel(lab[i] + '15' + '\n' + u[i])
            # start, end = ax[i].get_xlim()
            # ax[i].set_xticks(np.linspace(start, end, 1))
            # ax[i].legend(loc = 'center left')
            ax[i].xaxis.set_minor_locator(MultipleLocator(.1))
            ax[i].yaxis.set_minor_locator(MultipleLocator(minor[i]))
            ax[i].axvline(x=mc_t_pass-40., c='r', label = 'MW passage time')
            ax[i].axvline((pellet[n]['plasma_t'] - 40.), c='r', ls='--', lw=1, label = 'In-plasma time')

plt.xlabel("Time[s]")
#
ax01 = ax[0].twinx()
ax01.plot(pulse[n]['Bt_t']-40.,pulse[n]['Bt'],c='y')
ax01.set_ylabel('Bt' + '\n' + '[T]')
ax01.yaxis.set_minor_locator(MultipleLocator(.2))

# ax[0].yaxis.tick_right()
# ax[0].yaxis.set_label_position("right")
# ax01.legend(loc = 'center left')

ax[2].yaxis.set_minor_locator(MultipleLocator(.1))

# ax[-1].set_ylabel('Br n=1' + '\n'+ '[mT]')

axl1 = ax[-1].twinx()
axl1.plot(pulse[n]['my_phi_n1_t']-40.,pulse[n]['my_phi_n1'],c='y')
axl1.set_ylabel('phi_n1' + '\n'+ '[deg]')
axl1.set_ylim(-180,180)
axl1.set_yticks(np.linspace(-180, 180,4))
axl1.yaxis.set_minor_locator(MultipleLocator(30))
# axl1.legend(loc = 'center left')

plt.xlim(14,18)


e = time.time()
print('Computational Plotting time  is {:.3f}s'.format(e-s))

plt.show()



