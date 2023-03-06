import numpy as np
import time
from ppf import ppfdata
from pulses_datacollection import pseq_stef
from pulses_analysis import r_asy, pulse
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

start = time.time()


tstop = 0.015 #in seconds
s = int(2) #Harmonic number for ECE neasurements
B_tol = 0.01
qtol = 1e-6

#compute the radius at which we have q=2 +- 5%
for p in pseq_stef:
    totlenq  = len(ppfdata(p, 'EFIT', 'QMAG')[0])
    rowsq = int(totlenq)/33
    tq = ppfdata(p, 'EFIT', 'QMAG')[2]
    t0 =  tq[0]  # seconds
    tq = list(tq)

    tmax = r_asy[p]['tmax']
    l296 = 0
    while ppfdata(p, 'EFIT', 'QMAG')[1][l296] <2.96:
        l296+=1
    # now I retrieve the radius for (tmax - t0) onward (t0 is a time parameter which I can change
    # first of all I need to find the index at
    i = 0
    while tq[i] < (tmax - t0):
        i += 1
    j = i
    while tq[j] < (tmax - tstop):
        j += 1

    #Here I slice the matrix data into a smaller matrix qhich contains the values fron t0 to max in time and r>2.96 in space
    q_aux = ppfdata(p, 'EFIT', 'QMAG', reshape=(rowsq, 33))[0]
    q = q_aux[i:j,l296:]
    pulse[p]['q'] = q
    R_q = ppfdata(p, 'EFIT', 'QMAG')[1][l296:]
    pulse[p]['R_q'] = R_q
    colsq = len(R_q)
    twos = 2. * np.ones(len(R_q))
    R_q2 = []
    for k in range(len(q)):
        mq = min(abs(q[k]-2.))
        condition = np.allclose(q[k],twos)
        if condition is True:
            R_q2.append(0)
        else:
            for l in range(colsq):
                if (abs(q[k,l]-2.)-mq) < qtol:
                    R_q2.append(R_q[l])

    pulse[p]['R_q2'] = R_q2
    pulse[p]['R_q2_t'] = tq[i:j]
    #I want the index of B_toroidal at 2.96m
    Bt_296 = []
    totlenB = len(ppfdata(p, 'EFIT', 'BTAX')[0])
    rowsB = int(totlenB/33)
    Bt = ppfdata(p, 'EFIT', 'BTAX', reshape = (rowsB,33))[0]
    R_B = ppfdata(p, 'EFIT', 'BTAX')[1]

    m = min(abs(R_B-2.96))
    k296=0
    while abs((R_B[k296]-2.96) - m) > B_tol:
         k296+=1

    # I create a new radius vector which is the original one sliced from the (tmax - t0) to (tmax - tstop)
    # so from 1 s to 50 ms before the peak in  bolometry radiation, ie safetly before the disruption.
    Bt_296_dt = abs(Bt[i:j, k296])
    R_B_dt = R_B[k296]
    #Now I can calculate the frequency of the channel

    coef = s*28
    Fch = []
    for a in range(len(R_q2)):
        f = coef*R_B_dt*(Bt_296_dt[a])/R_q2[a]
        Fch.append(f)

    f_kk3 = ppfdata(p, 'KK3', 'GEN', reshape = (20,96))[0][15,:]
    Nch=[]

    for fr in Fch:
        c=0
        fmin = min(abs(f_kk3-fr))
        for channel, freq in enumerate(abs(f_kk3-fr)):
            if abs(freq - fmin) <1e-6:
                c = channel + 1
        Nch.append(c)

    pulse[p]['q2_chs'] = Nch

    # ch = (sum(Nch)/len(Nch))
    # r = ch - int(ch)
    # if r>0.5:
    #     ch = int(ch) + 1
    # else:
    #     ch  = int(ch)

    # pulse[p]['q2_kk3_ch'] = ch


#
# for p in pseq_stef:
#     print(p, pulse[p]['q2_kk3_ch'])
#
n = np.random.choice(pseq_stef)
# n=98071
#
# f, ax = plt.subplots(nrows=2, figsize = (9,12), sharex= True)
# f.suptitle('Pulse {:.0f}'.format(n))
#
# ax[0].plot(pulse[n]['R_q2_t'],pulse[n]['R_q2'], 'b-', ms=1)
# ax[0].set_ylabel('Major radius position of Q=2 surface [m]')
# ax[1].plot(pulse[n]['R_q2_t'],pulse[n]['q2_chs'], 'b-', ms=1)
# ax[1].set_ylabel('KK3 channel of Q=2 surface')
# # ax[1].axhline(pulse[n]['q2_kk3_ch'],c='r', ls ='--', label = 'Rounded up average channel = {}'.format(pulse[n]['q2_kk3_ch']))
# plt.xlabel('time [s]')
# plt.minorticks_on()
# # plt.legend()
# plt.show()
#
# for p in pseq_stef:
#     pulse[p]['Te_q2'] = ppfdata(p,'KK3', 'TE'+str(pulse[p]['q2_kk3_ch'])) [0] / (1e3)  # in keV
#     pulse[p]['Te_q2_t'] = ppfdata(p, 'KK3', 'TE' + str(pulse[p]['q2_kk3_ch'])+ '_t')[2]  # time vector
#
# for i in range(len(pulse[n]['q'])):
#     plt.plot(pulse[n]['R_q'],pulse[n]['q'][i],lw = 0.5)
# # plt.plot(pulse[n]['Te_q2_t'],pulse[n]['Te_q2'])
# plt.show()
end =time.time()
print('Computational q = 2 time is {:.3f}s'.format(end-start))