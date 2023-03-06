
from ppf import ppfdata
from pulses_analysis import r_asy
from pulses_pellets import pellet #for plotting the plasma arrival times
from q2 import pulse
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import sys
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

sys.path.append('KL8_files')
os.chdir('KL8_files')
plt.rcParams.update({'font.size':18})

shots = [98071,96103]

for n in shots:
    try:
        # q = 2 SURFACE TEMPERATURE RIGHT BEFORE THE DISRUPTION
        # RAW MICROWAVE SIGNAL
        pulse[n]['MCRW'] = ppfdata(n, 'MCRW', 'RAWS', uid='md7913')[0]  # signal in volts
        pulse[n]['MCRW_t'] = ppfdata(n, 'MCRW', 'RAWS', uid='md7913')[2]  # time vector
        # KL8-E8WA CAMERA INTEGRATED DATA (FROM JUVIL)
        pulse[n]['KL8-E8WA'] = np.genfromtxt('KL8-E8WA-' + str(n) + '-sum-1.txt', usecols=(1)) / 1e5
        pulse[n]['KL8-E8WA_t'] = np.genfromtxt('KL8-E8WA-' + str(n) + '-sum-1.txt', usecols=(0))
    except:
        OSError


q = ['Ip','Te92','powv', 'powh']
# q = ['Ip','Te92','Te_q2','nef3', 'nei3']
q1 = ['MCRW','nef3','loca']
q2 = ['powv', 'powh', 'KL8-E8WA']
q3 = ['MCRW','KL8-E8WA']
q4 = ['Ip','Te92','nef3']

lab = ['Ip', 'Te core', 'powV','powH']

u = ['[MA]','[keV]','[GW]', '[GW]']
# u = ['[MA]','[keV]','[keV]', '[10e20 m^-2]','[10e20 m^-2]']
u1= ['[V]','[10e20 m^-2]','[mT]']
u2 = ['[GW]', '[GW]', '[1e5 counts]']
u3 = ['[V]','[1e5 counts]']
u4= ['[MA]','[keV]','[10e20 m^-2]']

#NB : I HAVE MULTIPLIED EVERY TIME VECTOR BY 1000 TO DISPLAY IN MILLISECONDS
minor = [.5,.2,.1,.1,.1]

f, ax = plt.subplots(len(q), sharex=True, figsize= (16,9), facecolor = 'white')
f.suptitle(str(shots[0]))

col = ['r', 'b']
ypos = [1.03, 0.03]
for i, quant in enumerate(q):
    for j,n in enumerate(shots):
        ax[i].set_ylabel(lab[i] + '\n' + u[i])
        ax[i].set_xlim(0,100)
        #some pellets are broken, so I need to annotate two plasma times
        if type(pellet[n]['plasma_t']) is not list: # not broken
            t0 = pellet[n]['MW_t']
            ax[i].yaxis.set_minor_locator(MultipleLocator(minor[i]))
            ax[i].xaxis.set_minor_locator(MultipleLocator(1))
            # print(n,quant)
            ax[i].plot((pulse[n][quant + '_t']-t0)*1e3, pulse[n][quant],c = col[j],ls = '-', label='# '+ str(n), lw=1)
            ax[i].axvline((pellet[n]['plasma_t']- t0)*1e3, c=col[j], ls='--', lw = 1)#,label= str(n) + ' Plumes time to plasma')
            ax[0].annotate('+{:.1f}ms'.format(float((pellet[n]['plasma_t']- t0)*1e3)), xy = ((pellet[n]['plasma_t']- pellet[n]['MW_t'])*1e3,ypos[j]),
                            xycoords = ('data','axes fraction'))
            ax[i].axvline(25.3, c='b', ls='--', lw=1)
            ax[0].annotate('+25.3ms', xy = (25.3,1.03),
                            xycoords = ('data','axes fraction'))
        else:
            print(12)
            t0 = pellet[n]['MW_t']
            ax[i].plot((pulse[n][quant + '_t']-t0)*1e3, pulse[n][quant],c = col[j],ls='-', label='# '+ str(n), lw=1)
            ax[i].yaxis.set_minor_locator(MultipleLocator(minor[i]))
            for shot,t in enumerate((pellet[n]['plasma_t']- t0)*1e3):
                print(t)
                ax[i].axvline(t, c = col[j],ls='--', lw=1)#, label = str(n) + ' Plumes time to plasma' )
                ax[0].annotate('+{:.1f}s'.format(float(t)), xy=(t,ypos[j]),xycoords=('data', 'axes fraction'))
            # ax[i].xaxis.set_minor_locator(MultipleLocator(.1))
            # ax[i].xaxis.set_major_locator(MultipleLocator(5))

ax[0].legend(loc='upper right', framealpha=1)
plt.xlabel("Time after first peak in MCRW camera signal [ms]")
plt.show()