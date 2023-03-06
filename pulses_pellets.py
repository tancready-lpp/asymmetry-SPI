# coding=utf-8

from pulses_datacollection import pseq_stef
import time
import numpy as np
import ppf
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from getdat import getdat8

def sign(x):
    M = max(x)
    m = min(x)
    t = 1
    if m<-t:
        s=-1
    elif M>t:
        s=1
    else:
        s = 0
    return s



# pellet_type = ['100%Ne', '100%Ne', '100%Ne', '100%Ne',  '100%Ne', '100%Ne', '100%Ne', '100%Ne', '100%Ne',
#                '100%Ne', '100%Ne', '100%Ne', '100%Ne', '100%Ne', '10%Ne', '10%Ne'] #'100%Ne',

t_to_plasma = np.genfromtxt('plasma_phi1.txt', usecols=1)

#pulse = {
          #efcc15 #efcc37 #pellet type #Barrel #MW Cavity Time (s)  #Time to Plasma (s) #Pellet Speed (m/s) #N pieces

start = time.time()

barrel_lenght =  4.852 # from half MW camera to SOL # 4.715m from vessel top
user = 'md7913'

pellet = dict()
for i, p in enumerate(pseq_stef):
    pellet[p] = dict()
    pellet[p]['efcc15'] = sign(getdat8("C2E-EFCC<15", p)[0]  / 1e3)
    pellet[p]['efcc37'] = sign(getdat8("C2E-EFCC<37", p)[0]  / 1e3)
    # pellet[p]['type'] = pellet_type[i]
    pellet[p]['barrel'] = 'B'

    # searching the time at which we have a relative minimum in MW signal
    V_mw = abs(ppf.ppfdata(p, 'MCRW', 'FIT', uid=user)[0]) # if the pellet is broken we'll have multiple values
    peaks_ind = argrelmax(V_mw)
    pellet[p]['MW_t'] = ppf.ppfdata(p, 'MCRW', 'FIT', uid=user)[2][peaks_ind]
    pellet[p]['pieces'] = (ppf.ppfdata(p, 'MCRW', 'NPEK', uid=user)[0])
    pellet[p]['plasma_t'] = t_to_plasma[i]

    #if pellet is not broken after MW
    if type((pellet[p]['plasma_t'])) is not list:
    # This if statement is needed since some pulses have no camera data for t_to_plasma so I impose it to be the one computed
    # in the ppf library (NB: could be very distant from the real value)
        if t_to_plasma[i] <1:
            pellet[p]['comment'] = 'estimated'
            pellet[p]['velocity'] = ppf.ppfdata(p, 'MCRW', 'VLCT', uid = user)[0][0]
            pellet[p]['plasma_t'] = barrel_lenght/pellet[p]['velocity'] + float(pellet[p]['MW_t'])
        else:
            pellet[p]['comment'] = 'good'
            pellet[p]['velocity'] = barrel_lenght / (pellet[p]['plasma_t'] - (pellet[p]['MW_t']))
    #If pellet is broken after MW
    else:
        pellet[p]['velocity'] = barrel_lenght / (pellet[p]['plasma_t'] - (pellet[p]['MW_t']))
        pellet[p]['pieces'] = len((pellet[p]['plasma_t']))
        pellet[p]['comment'] = 'maybe broken after MW'


end = time.time()
print('Computational Pellet info time = {:.3f}s'.format(end-start))