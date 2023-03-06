from pulses_datacollection import pulse, pseq_stef
from pulses_analysis import raf
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import trapezoid
from scipy.stats import chisquare
#in this script I want to fit the rafmax_VS_phi on a cosine distribution and try to calculate and model the TPF

k = (180./np.pi) #from radiants to degrees

def p_profile(phi, dp,dphi, l,d):
    n = np.exp(-(phi/l)**2.)
    p = 1.+dp*np.cos(d-dphi-phi)
    pi = np.pi
    x = np.linspace(-pi, +pi, 360)
    avg1 = trapezoid(p * n, x, dx=1e-6) / (2 * pi)
    n0 = 1. / avg1
    n *= n0
    avg = trapezoid(p * n, x, dx=1e-6) / (2 * pi)
    P=n*p*avg
    return P


def prad_func(phi, dp, dphi,l):
    spi = 11.25/k
    phi_n1_formula = phi + spi
    # vloc = 90./k - spi
    # hloc = -135./k - spi
    vloc = 90./k
    hloc = -135./k
    pv = 1.+ dp*np.cos(phi_n1_formula - dphi - vloc)
    ph = 1.+ dp*np.cos(phi_n1_formula - dphi - hloc)
    nv = np.exp(-((vloc-spi)/l)**2.)
    nh = np.exp(-((hloc-spi)/l)**2.)
    Pv = pv*nv
    Ph = ph*nh
    raf = (Pv-Ph)/(Pv+Ph)
    return raf


def tpf_calc(phi_n1, dp, dphi,l):
    pi = np.pi
    a = -pi
    b = pi
    angle = np.linspace(a,b,360)

    p = 1. + dp*np.cos(phi_n1 + dphi - angle)
    n = np.exp(-(angle/l)**2.)
    x = np.linspace(-pi, +pi, 360)
    #n0 computation from integral
    avg = trapezoid(p*n, x,dx=1e-6)/(2*pi)
    n0 = 1./avg
    n*=n0
    tpf = max(p*n)
    return tpf

#--------MAIN--------#
plt.rcParams.update({'font.size':15})
#---RAF---#
xdata = [x/k for x in raf['phi_n1']]  #in radiants



ydata = raf['rafmax']
yerr = raf['err_rafmax']


s = [y*0.2 for y in ydata]


#fitting the raf_func
popt, pcov = curve_fit(prad_func, xdata,ydata, sigma = s)
perr = np.sqrt(np.diag(pcov))
# print(popt[0],perr[0], popt[1]*k, perr[1]*k, popt[2]*k, perr[2]*k)
ydataplus = list()
ydataminus = list()
for i in range(len(ydata)):
    ydataplus.append(ydata[i] + yerr[i])
    ydataminus.append(ydata[i] - yerr[i])

splus = [y*0.2 for y in ydataplus]
smin = [y*0.2 for y in ydataminus]

pplus,covplus = curve_fit(prad_func,xdata, ydataplus, sigma = splus)
pmin,covmin = curve_fit(prad_func,xdata, ydataminus, sigma = smin)

if popt[0]<0:
    popt[0]*=-1
    popt[1]+=np.pi


print('\nraf_func fit STATISTICAL\ndp = {:.2f} pm {:.2f}\ndphi = ({:.1f} pm {:.1f})deg\nlambda = ({:.1f} pm {:.1f}) deg'.format(popt[0],perr[0], (popt[1])*k,perr[1]*k, popt[2]*k, perr[2]*k))

#
# print('OHMIC\nPOPT\ndp = {:.2f}\ndphi = {:.2f}\nl={:.2f}'.format(popt[0],popt[1]*k-360,popt[2]*k))
# print('\nPPLUS\ndp = {:.2f}\ndphi = {:.2f}\nl={:.2f}'.format(pplus[0],pplus[1]*k,pplus[2]*k))
# print('\nPMIN\ndp = {:.2f}\ndphi = {:.2f}\nl={:.2f}\n'.format(-pmin[0],pmin[1]*k-180,pmin[2]*k))

# print('\nraf_func fit DATA +- ERR\ndp = {:.2f} pm {:.2f}\ndphi = ({:.1f} pm {:.1f})deg\nlambda = ({:.1f} pm {:.1f}) deg'.format(popt[0],perr[0], (popt[1])*k,perr[1]*k, popt[2]*k, perr[2]*k))



phi = np.linspace(-90,270,360)
# phi = np.linspace(-180,180,360)
raf_fit = prad_func(phi/k, popt[0],popt[1],popt[2])
raf_fit_errpos = prad_func(phi/k, pplus[0], pplus[1], pplus[2])
raf_fit_errneg = prad_func(phi/k, pmin[0], pmin[1], pmin[2])

raf_dp_p = prad_func(phi/k,popt[0]+perr[0], popt[1],popt[2])
raf_dp_m = prad_func(phi/k,popt[0]-perr[0], popt[1],popt[2])

raf_dphi_p = prad_func(phi/k,popt[0], popt[1]+perr[1],popt[2])
raf_dphi_m = prad_func(phi/k,popt[0], popt[1]-perr[1],popt[2])

raf_l_p = prad_func(phi/k,popt[0], popt[1],popt[2]+perr[2])
raf_l_m = prad_func(phi/k,popt[0], popt[1],popt[2]-perr[2])

# --- CHI2 SQUARE TEST --- #
# raf_chi2 =[]
# for x in xdata:
#     ind = 0
#     while abs(phi[ind]/k - x)>1e-2:
#         ind+=1
#     raf_chi2.append(raf_fit[ind])
# chi2 = chisquare(ydata, raf_chi2)
# print('Chi2 = ', chi2[0])
# print('p-value = ', chi2[1])


l_s = 135./k
dphi_s = -20./k
dp_s = 0.22
fit_s = prad_func(phi/k,dp_s,dphi_s,l_s)

#---TPF---#
phi_rad = phi/k
tpf = []
tpf_s = []
tpf_errpos =[]
tpf_errneg = []

tpf_dp_p = []
tpf_dp_m = []
tpf_dphi_p = []
tpf_dphi_m = []
tpf_l_p = []
tpf_l_m = []
tpf_ideal = []
tpf_non = []

for phi_n1 in phi_rad:
    tpf.append(tpf_calc(phi_n1, popt[0], popt[1],popt[2]))
    tpf_errpos.append(tpf_calc(phi_n1, pplus[0], pplus[1], pplus[2]))
    tpf_errneg.append(tpf_calc(phi_n1, pmin[0], pmin[1], pmin[2]))
    tpf_s.append(tpf_calc(phi_n1, dp_s, dphi_s, l_s))
    tpf_dp_p.append(tpf_calc(phi_n1,popt[0]+perr[0],popt[1],popt[2]))
    tpf_dp_m.append(tpf_calc(phi_n1,popt[0] - perr[0], popt[1], popt[2]))
    tpf_dphi_p.append(tpf_calc(phi_n1,popt[0] , popt[1]+ perr[1], popt[2]))
    tpf_dphi_m.append(tpf_calc(phi_n1,popt[0], popt[1] - perr[1], popt[2]))
    tpf_l_p.append(tpf_calc(phi_n1,popt[0], popt[1] , popt[2]+ perr[2]))
    tpf_l_m.append(tpf_calc(phi_n1,popt[0], popt[1], popt[2] - perr[2]))
    tpf_ideal.append(tpf_calc(phi_n1, 0.1, 0, 175/k))
    tpf_non.append(tpf_calc(phi_n1, 0.45, 0, 100/k))

dp_aux = np.linspace(0,0.5,50)
l_aux = np.linspace(90/k,180/k,90)

m_tpf = max(tpf)
m_tpf_dp_p = max(tpf_dp_p)
m_tpf_dp_m = max(tpf_dp_m)
m_tpf_l_p = max(tpf_l_p)
m_tpf_l_m = max(tpf_l_m)


i=0
while abs(tpf[i]-m_tpf)>1e-4:
    i+=1

loc_max = phi_rad[i]*k
print('location Max tpf = {:.2f} deg'.format(loc_max))

print('TPF max = {:.3f}'.format(m_tpf))
# print('TPF max dp min = {:.3f}'.format(m_tpf_dp_m))
# print('TPF maxdp plus = {:.3f}'.format(m_tpf))
# print('TPF max l min = {:.3f}'.format(m_tpf_l_m))
# print('TPF max l plus= {:.3f}'.format(m_tpf_l_p))

# --------2d map ---------- #

# TPF = []
# max_tpf = np.zeros((len(dp_aux),len(l_aux)),float)
#
# dphi = 0
#
# for i,dp in enumerate(dp_aux):
#     for j,l in enumerate(l_aux):
#         for phi_n1 in phi_rad:
#             TPF.append(tpf_calc(phi_n1,dp,dphi,l))
#         max_tpf[i,j] = max(TPF)
#         TPF = []
#
# max_tpf = np.transpose(max_tpf)
#
# print('TPF max ideal = {:.3f}'.format(max(tpf_ideal)))
# print('TPF max non-ideal = {:.3f}'.format(max(tpf_non)))
#
#
# plt.title('Max TPF')
# plt.contourf(dp_aux,l_aux*k,max_tpf,levels = 50, cmap = 'jet')
# plt.colorbar()
# cs = plt.contour(dp_aux,l_aux*k,max_tpf,levels = 50, colors = 'black')
# plt.clabel(cs,cs.levels[::3], inline=1, fontsize = 12)
# plt.plot(0.18, 166, marker ='p',c='k', ms = 10,ls='None', label = 'H-mode set')
# plt.plot(0.36, 129, marker = '*', c='k', ms = 10,ls='None', label = 'Ohmic set')
# plt.plot(0.1, 175, marker = 'o', c='k', ms = 10,ls='None', label = 'Ideal case')
# plt.plot(0.45, 100, marker = 's', c='k', ms = 10,ls='None', label = 'Not-ideal case')
# plt.legend(loc=0)
# plt.minorticks_on()
#
#
# # print(max_tpf)
# # plt.plot(dp_aux,m_tpf_dp,'bo')
# plt.xlabel('$\Delta p$')
# # plt.plot(l_aux*k,m_tpf_l,'bo')
# # plt.xlabel('lambda')
# plt.ylabel('$\lambda_{\phi}$ [deg]')
# plt.show()
# # #

# ----power profile
#
# d = 30/k
# #
# phi_aux = np.linspace(-180,180,360)/k
# P_h1= p_profile(phi_aux,0.19,-13/k,163.2/k,d)
# P_ohm= p_profile(phi_aux,0.34,-18/k,126.7/k,d)
# ideal = p_profile(phi_aux,0.1,0, 175/k,d)
# worst = p_profile(phi_aux,0.45,0, 100/k,d)
# # #
# # #
# # # l=163.2/k
# # #
# # # g = np.exp(-((phi_aux)/l)**2)
# # #
# plt.plot(phi_aux*k,ideal, 'k--', label = 'Ideal')
# plt.plot(phi_aux*k,worst, 'g--', label = 'Not-Ideal')
# plt.plot(phi_aux*k,P_h1, 'r-', label = '$H-mode$')
# plt.plot(phi_aux*k,P_ohm, 'b-', label = 'Ohmic')
# plt.xlabel('Toroidal angle respct to SPI position [deg]')
# plt.ylabel('$P_{rad}$ [GW]')
# plt.xticks(np.arange(-180,180,30))
# # plt.xlim
# plt.title('Power profile')
# plt.minorticks_on()
#
# plt.legend()
#
# plt.show()




s_phi,s_tpf,s_raf = np.genfromtxt('paper_stefan_data.txt', usecols = (0,1,2), unpack=True)


# #----- PLOTTING -----#
plt.rcParams.update({'font.size':18})
plt.figure(figsize=(12,7), facecolor='white')
plt.xticks(np.arange(-90,360,45))
plt.xlim(-90,270)
plt.xlabel('Phi_n=1 locked mode O-point phase [deg] respect to SPI angular position')


# #--- RAF PLOTTING --- #

# plt.fill_between(phi, raf_fit_errpos, raf_fit_errneg,color= 'cornflowerblue')
plt.fill_between(phi, raf_dp_p, raf_dp_m,color= 'orange')
plt.fill_between(phi, raf_l_p, raf_l_m,color= 'orange')
plt.fill_between(phi, raf_dphi_p, raf_dphi_m,color= 'orange')
#
# # plt.plot(phi,raf_fit_errpos,ls = 'dotted', c='b')
plt.plot(raf['phi_n1'], raf['rafmax'],'bo', markersize= 10, label = 'dataset')
# # plt.plot(raf['my_phi_n1_tmax'], raf['rafmax'],'ro', markersize= 10, label = 'my dataset')
# # plt.plot(raf['my_phi_n1_ppf_tmax'], raf['rafmax'],'go', markersize= 10, label = 'ppf dataset')
plt.plot(phi,raf_fit, 'b-', label = 'Fit of RAF')
# # plt.plot(phi,raf_fit_errneg,ls = 'dotted', c='b', label = 'Fit of RAF $\pm \Delta$RAF')
plt.errorbar(raf['phi_n1'], raf['rafmax'],xerr = raf['dphi_tmax'],yerr=raf['err_rafmax'],fmt = 'none', c='b', capsize=3)
plt.axvline(x=0., ls='--', c='y', label = 'SPI location')

# plt.plot(phi,raf_l_p,ls = 'dotted', c='r', label = 'Raf fit pm covariance errors on lambda')
# plt.plot(phi,raf_l_m,ls = 'dotted', c='r')
# plt.plot(phi,raf_dp_p,ls = 'dotted', c='g', label = 'Raf fit pm covariance errors on dp')
# plt.plot(phi,raf_dp_m,ls = 'dotted', c='g')
# plt.plot(phi,raf_dphi_p,ls = 'dotted', c='k', label = 'Raf fit pm covariance errors on dphi')
# plt.plot(phi,raf_dphi_m,ls = 'dotted', c='k')
# plt.plot(phi,raf_fit_errpos, ls = 'dotted', c='b')
# plt.plot(phi,raf_fit_errneg, ls = 'dotted', c='b', label = 'Fit + error range')
# plt.plot(phi,fit_s,'r--', label = 'Fitted RAF curve with Stefan fit parameters')
# plt.plot(s_phi,s_raf,'g-', label = 'Paper RAF Stefan')
#
plt.minorticks_on()

# plt.text(-89,0.1, 'Chi2 = {:.3f}\np-value = {:.3f}'.format(chi2[0], chi2[1]), weight = 'bold')
#
# txt = 'My FIT (H-mode)\n$\Delta$p = {:.2f} $\pm$ {:.2f}\n$\Delta \phi_1$ = {:.1f}° $\pm$ {:.1f}°\n$\lambda_n$ =  {:.1f}° $\pm$ {:.1f}°'.format(popt[0],perr[0],(popt[1])*k, perr[1]*k,popt[2]*k, perr[2]*k)
# plt.text(0,-2, s=txt)
#
# txt_s = 'Stefan FIT (H-mode)\n$\Delta$p_s = {:.3f}\n$\Delta \phi_1s$ = {:.1f}°\n$\lambda_ns$ =  {:.1f}°'.format(dp_s,dphi_s*k,l_s*k)
# plt.text(50,0.2, s=txt_s)

plt.ylabel('Radiation Asymmetry Fraction at tmax')
plt.legend(loc='best', fontsize = 9, markerscale = 0.75)

# for i, p in enumerate(pseq_stef):
#     plt.annotate('{:.0f}'.format(p), (raf['phi_n1'][i], raf['rafmax'][i]), textcoords = 'offset points', xytext =(5,5))
plt.show()

# # #--- TPF PLOTTING --- #
plt.figure(figsize=(12,7), facecolor='white')
plt.xticks(np.arange(-90,360,45))
plt.xlim(-90,270)
plt.xlabel('Phi_n=1 locked mode O-point phase [deg] respect to SPI angular position')
#
# plt.suptitle('Toroidal Peacking Factor VS Phase of n=1 locked mode for H-mode plasma')
# plt.fill_between(phi, tpf_errpos, tpf_errneg, color='cornflowerblue')
plt.fill_between(phi, tpf_dp_p, tpf_dp_m, color='orange')
plt.fill_between(phi, tpf_dphi_p, tpf_dphi_m, color='orange')
plt.fill_between(phi, tpf_l_p, tpf_l_m, color='orange')
# plt.plot(phi,tpf_errpos, c='b', ls = 'dotted', label = 'TPF limits')
# plt.plot(phi,tpf_errneg, c='b', ls = 'dotted')
plt.plot(phi, tpf, 'b-', label = 'TPF')
# # plt.plot(phi, tpf_s, 'r--', label = 'TPF Stefan')
# # plt.plot(s_phi,s_tpf,'g-', label = 'Paper stefan TPF')
plt.axvline(x=0., ls='--', c='y', label = 'SPI location')
plt.axvline(x=loc_max, ls='dotted', c='r')
plt.axhline(y=max(tpf), c = 'r',ls='dotted')
plt.annotate('TPF max = {:.2f} at {:.2f} deg'.format(m_tpf, loc_max), (125,m_tpf),textcoords = 'offset points', xytext =(2,2))
plt.ylabel('Toroidal Peaking Factor')
plt.legend(loc='lower left')
plt.minorticks_on()
# #
#
#
# plt.show()