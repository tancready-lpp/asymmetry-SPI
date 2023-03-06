# Created 12:02 Wed 24 Aug 2022

# This code aims to do store the data needed for STEP1 of my work plan in a nesterd dictionary. 
# I will write all the comments as I code.

# The following libraries allow me to enter and retreive JET data, processed or not

import numpy as np
from ppf import ppfdata  # this if the library for Processed Pulse Data, and relative useful methods (ppf)
from getdat import getdat8  # this if the library for Jet Pulse Data and relative useful method(jpf)
import time
import matplotlib.pyplot as plt



def phierror(sx01,sx14,sy01,sy14,s101,s114,s501,s514, l):
    data = [sx01, sx14, sy01, sy14, s101, s114, s501, s514]
    data = np.array(data)
    coef = [0.0524,0.0534,-0.0536,-0.0526,0.053,0.0537,-0.053,-0.0523]
    sin = np.zeros(len(sx01))
    cos = np.zeros(len(sx01))
    for k in range(len(sin)):
        for i in range(3):
            sin[k] += (data[i][k] * coef[i])*0.25
        for i in range(4,8):
            cos[k] += (data[i][k] * coef[i])*0.25

    phi = np.zeros(l)
    for k in range(l):
        if cos[k]!=0:
            phi[k] = np.arctan(sin[k]/cos[k])

    ddata = np.zeros((len(data),len(sin)))
    dsin = np.zeros(len(sin))
    dcos = np.zeros(len(sin))
    for k in range(len(sin)):
        for i in range(len(data)):
            ddata[i][k] = max(2e-3,0.015*data[i][k])

    for k in range(len(sin)):
        for i in range(3):
            dsin[k] += ddata[i][k]*abs(coef[i])/4

    for k in range(len(sin)):
        for i in range(4,8):
            dcos[k] += ddata[i][k]*abs(coef[i])/4

    dphi = np.zeros(l)
    amp = sin**2 + cos**2
    dphi = abs(cos/amp)*dsin + abs(sin/amp)*dcos

    return  phi,dphi

def phi_err(x,y):
    dx = 0.015*abs(x)
    dy = 0.015*abs(y)
    amp2 = x**2+y**2
    dphi = np.zeros(len(amp2))
    phi = np.zeros(len(amp2))
    # dphi = abs(x)*dy + abs(y)*dx
    for k in range(len(amp2)):
        if amp2[k] != 0 and y[k]!=0:
            phi[k] = np.arctan(x[k]/y[k])
            dphi[k] = abs(y[k]/amp2[k])*dx[k] + abs(x[k]/amp2[k])*dy[k]
        else:
            dphi[k] = 0

    return phi, dphi


# ------- MAIN -------- #

start = time.time()


pseq_stef, phi_n1 = np.genfromtxt('plasma_phi1.txt', usecols = (0,2), unpack=True) #for OHMIC MODE

# pseq_stef, phi_n1 = np.genfromtxt('h-mode.txt', usecols = (0,1), unpack=True) #for H-MODE


# I store the variables in a nested dictiorary; the main one will be the pulse,
# and the nested one will have the names of the variables as keys and lists/arrays as values.
# Some of the varialbles will be retreived form jpf while others from ppf.

# An example of the data nested dictionary is:
# pulse[pulse number] = {
    # Ip: [nparray] = plasma current - from jpf
    # Ip_t: [nparray] = time variable for plasma current
    # nei3: [nparray, len = 701] = integrated electron density from channel 3 vertical through the core - from ppf
    # nei3_t: [nparray, len = 701] = time variable for integrated electron density from channel 3 vertical
    # through the core - from ppf
    # efcc15: [nparray, len = 34500] = efcc current octant 1-5 - from jpf
    # efcc15_t: [nparray, len = 34500]  =time vector for efcc current octant 1-5 - from jpf
    # efcc37: [nparray, len = 34500] = efcc current octant 3-7 - from jpf
    # efcc37_t: [nparray, len = 34500] = time vector for efcc current octant 3-7 - from jpf
    # loca: [nparray, len = 12068] = locked mode (n=1) amplitude - from jpf
    # loca_t: [nparray, len = 12068] = time vector for locked mode (n=1) amplitude - from jpf
    # phi_n1: [list I think, len = 17]  = locked mode n=1 phase in SPI frame of reference [degrees] - from excel
    # powv: [nparray, len = 106736] = Vertical Bolometer Power data during disruption - from another UID! (ppf)
    # powv_t: [nparray,len = 106736] = Time vector for Vertical Bolometer Power data during disruption
    # powh:  [nparray,len = 106736] = Horizontal Bolometer Power data during disruption -from another UID! (ppf)
    # powh_t: [nparray,len = 106736] = Time vector for Horizontal Bolometer Power data during disruption
    #   }

pulse = dict()

for n in pseq_stef:
    pulse[n] = dict()

for i, n in enumerate(pseq_stef):
    # PLASMA CURRENT
    pulse[n]['Ip'] = getdat8('C2-IPLA', n)[0] * ((-1) / (1e6))  # making the current signal positive and in MA
    pulse[n]['Ip_t'] = getdat8('C2-IPLA', n)[1]  # time vector
    # TOROIDAL MAGNETIC FIELD
    pulse[n]['Bt'] = ppfdata(n, 'MAGN', 'BVAC')[0]*(-1) # Making the toroidal filed positive
    pulse[n]['Bt_t'] = ppfdata(n, 'MAGN', 'BVAC')[2]  # time vector
    # ELECTRON TEMPERATURE FROM ECE INTERFEROMETER CH 92 (CORE) and 32 (EDGE)
    pulse[n]['Te92'] = ppfdata(n,'KK3', 'TE92') [0] / (1e3)  # in keV
    pulse[n]['Te92_t'] = ppfdata(n, 'KK3', 'TE92')[2] #time vector
    # pulse[n]['Te32'] = ppfdata(n, 'KK3', 'TE32')[0] / (1e3)  # in keV
    # pulse[n]['Te32_t'] = ppfdata(n, 'KK3', 'TE32')[2]  # time vector
    # ELECTRON DENSITY ON MAGNETIC AXIS
    pulse[n]['nei3'] = ppfdata(n,'KG1V','LID3')[0] / (1e20)  # in units of 10^20
    pulse[n]['nei3_t'] = ppfdata(n, 'KG1V', 'LID3')[2]  # time vector
    # ELECTRON DENSITY FROM FARADAY EFFECT (channel 3)
    pulse[n]['nef3'] = getdat8('G4R-LID:003', n)[0] / (1e20)
    pulse[n]['nef3_t'] = getdat8('G4R-LID:003', n)[1]  # time vector
    # EFCC CURRENT FROM BOTH OCTANTS-COUPLE (1-5 & 3-7)
    pulse[n]['efcc15'] = getdat8("C2E-EFCC<15", n)[0]  / 1e3  # reporting in kA
    pulse[n]['efcc15_t'] = getdat8("C2E-EFCC<15", n)[1]  # time vector
    pulse[n]['efcc37'] = getdat8("C2E-EFCC<37", n)[0] / 1e3  # reporting in kA
    pulse[n]['efcc37_t'] = getdat8("C2E-EFCC<37", n)[1]  # time vector
    # LOCKED n=1 MODE AMPLITUDE AND PHASE
    pulse[n]['loca'] = getdat8("C2-LOCA", n)[0] / (1e-3)  # measured in mT
    pulse[n]['loca_t'] = getdat8("C2-LOCA", n)[1]
    cos = np.array(getdat8('C2-MHDG',n)[0])
    sin = np.array(getdat8('C2-MHDF',n)[0])

    # print(n,len(cos))
    tan = np.zeros(len(cos))
    for j in range(len(cos)):
        if cos[j] !=0:
            tan[j] = (sin[j])/(cos[j])
        else:
            tan[j]=0
    # pulse[n]['my_phi_n1'] = np.arctan(tan) * (180 / np.pi)
    pulse[n]['my_phi_n1_t'] = getdat8("C2-MHDF", n)[1]
    pulse[n]['my_phi_n1_ppf_t']  =  pulse[n]['my_phi_n1_t']


    pulse[n]['phi_n1'] = phi_n1[i]
    sx01 = getdat8('C2-SX01',n)[0]
    sx14 = getdat8('C2-SX14', n)[0]
    sy01 = getdat8('C2-SY01',n)[0]
    sy14 = getdat8('C2-SY14',n)[0]
    s101 = getdat8('C2-S101',n)[0]
    s114 = getdat8('C2-S114',n)[0]
    s501 = getdat8('C2-S501',n)[0]
    s514 = getdat8('C2-S514', n)[0]

    dphi = phierror(sx01,sx14,sy01,sy14,s101,s114,s501,s514, len(sin))[1]*(180/np.pi)
    # dphi = phi_err(sin,cos)[1] * (180/np.pi)


    # pulse[n]['my_phi_n1'] = phierror(sx01,sx14,sy01,sy14,s101,s114,s501,s514, len(sin))[0]*(180/np.pi) + 90 + 11.25
    pulse[n]['my_phi_n1_ppf'] = np.arctan(tan)*(180/np.pi) + 90
    pulse[n]['my_phi_n1'] = phi_err(sin,cos)[0] * (180/np.pi) + 47


    pulse[n]['dphi'] = dphi
    pulse[n]['dphi_t'] = getdat8('C2-MHDF',n)[1]
    # pulse[n]['my_phi_n1'] = np.arctan(getdat8('C2-MHDF',n)[0]/getdat8('C2-MHDG',n)[0]) * (180/np.pi)

    # POW VERTICAL
    pulse[n]['powv'] = ppfdata(n, 'RAD', 'POWV', uid='sjach')[0] / (1e9)  # measure in GW
    pulse[n]['powv_t'] = ppfdata(n, 'RAD', 'POWV', uid='sjach')[2] + 8e-4  # time delay of 0.8ms correction
    if len(pulse[n]['powv']) == 0:
        del pulse[n]['powv']
        del pulse[n]['powv_t']
        pulse[n]['powv'] = ppfdata(n, 'RAD', 'POWV', uid='mlehnen')[0] / (1e9)  # measure in GW
        pulse[n]['powv_t'] = ppfdata(n, 'RAD', 'POWV', uid='mlehnen')[2] + 8e-4   # time delay of 0.8ms correction
    # POW HORIZONTAL
    pulse[n]['powh'] = ppfdata(n, 'RAD', 'POWH', uid='sjach')[0] / (1e9)  # measure in GW
    pulse[n]['powh_t'] = ppfdata(n, 'RAD', 'POWH', uid='sjach')[2] + 8e-4   # time delay of 0.8ms correction
    if len(pulse[n]['powh']) == 0:
        del pulse[n]['powh']
        del pulse[n]['powh_t']
        pulse[n]['powh'] = ppfdata(n, 'RAD', 'POWH', uid='mlehnen')[0] / (1e9)  # measure in GW
        pulse[n]['powh_t'] = ppfdata(n, 'RAD', 'POWH', uid='mlehnen')[2] + 8e-4  # time delay of 0.8ms correction






# i=0
# for p in pseq_stef:
#     if np.nan in pulse[p]['Tcore']:
#         print(p)
#     else:
#         print(str(p)+' no')
#         i+=1
# print(i)
# print(len(pseq_stef))

end = time.time()
print('Computational Collection time  is {:.3f}s'.format(end-start))

