#XRD Powder Diffraction Analysis Program (for cubic crystal systems)
#Author = MV Mastwijk, 15-OCT-2021
#See below for the settings

from sys import path_importer_cache
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy.core.function_base import linspace
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths, argrelmin, detrend
from scipy.integrate import quad
import pathlib

from scipy.sparse import base

#settings

#Cu2O-A SETTINGS
filename = "cu2o_A-X_Y.txt" #note: file must be in the same folder as the program
cutoff_2xtheta = 12 #degrees
peak_prominence = 20 #higher peak prominence > small peaks are ignored

wavelength = 0.71073 #A^0
#

#import data
powder_diffraction = np.genfromtxt(pathlib.PurePath(pathlib.Path(__file__).parent, filename))
x2theta,intensity=[],[]
for i in range(len(powder_diffraction)):
    x2theta.append(powder_diffraction[i][0])
    intensity.append(powder_diffraction[i][1])

#setting up plot
fig, ax = plt.subplots(1, figsize=(8, 6), facecolor='whitesmoke')
ax.set_facecolor("gainsboro")
ax.grid(color='white', linestyle='-', linewidth=1)
fig.suptitle(filename, fontsize=15)
plt.xlabel("2θ (°)")
plt.ylabel("I (arbitrary units)")
#ax.plot(x2theta, intensity, color="orange")
#plt.axvline(x=cutoff_2xtheta, linestyle='--', color="grey")

#peak finding
peak_list, peak_data=find_peaks(intensity, prominence=peak_prominence, width=0, rel_height=1)
peak_width=peak_widths(intensity,peak_list,rel_height=0.6)[0]

#baseline
negintensity=[]
for i in range(len(intensity)):
    negintensity.append(-intensity[i])
baseline, baseline_data = find_peaks(negintensity)
baseline=np.append(baseline,[len(x2theta)-1])
basex2theta=[]
baseintensity=[]
for i in range(len(baseline)):
    #plt.plot(x2theta[baseline[i]], intensity[baseline[i]], ".", color="blue")
    basex2theta.append(x2theta[baseline[i]])
    baseintensity.append(intensity[baseline[i]])
baseline_curve=np.poly1d(np.polyfit(basex2theta,baseintensity,10))
#x=np.linspace(0,44,1000000)
#ax.plot(x, baseline_curve(x), color="red")
intensity_corrected=[]
for i in range(len(x2theta)):
    intensity_corrected.append(intensity[i]-baseline_curve(x2theta[i]))

x2theta_dummy=[]
intensity_dummy=[]
for i in range(len(x2theta)):
    if x2theta[i]>=cutoff_2xtheta:
        x2theta_dummy.append(x2theta[i])
        intensity_dummy.append(intensity_corrected[i])
ax.plot(x2theta_dummy, intensity_dummy, color="orange")

#finding initial parameters for fitting
peak_x,peak_y=[],[]
base_l,base_r=[],[]
n=0
for i in range(len(peak_list)):
    if x2theta[peak_list[i]]<cutoff_2xtheta:
        n+=1
    else:
        peak_x.append(x2theta[peak_list[i]])
        peak_y.append(intensity_corrected[peak_list[i]])
        peak_width[i]=int(round(peak_width[i]))
        base_l.append(x2theta[int(peak_list[i]-peak_width[i])])
        base_r.append(x2theta[int(peak_list[i]+peak_width[i])])
for i in range(n):
    peak_list=np.delete(peak_list, 0)
    peak_width=np.delete(peak_width, 0)

#fit peaks
def _1Voigt(x, ampG1, cen, sigmaG1, ampL1, widL1):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen)**2)/((2*sigmaG1)**2)))) +\
              ((ampL1*widL1**2/((x-cen)**2+widL1**2)) )

peak_position=[]
peak_parameters=[]
for i in range(len(peak_list)):
    intensity_range=[]
    x2theta_range=[]
    for j in range(len(x2theta)):
        if x2theta[j]>=base_l[i] and x2theta[j]<=base_r[i]:
            x2theta_range.append(x2theta[j])
            intensity_range.append(intensity_corrected[j])
    halfintensity = intensity_corrected[peak_list[i]]/2
    gauss_stddev = np.abs(x2theta[peak_list[i]]-base_l[i])
    lorentz_width = np.abs(x2theta[peak_list[i]]-base_l[i])*1000
    parameters, info = curve_fit(_1Voigt, x2theta_range, intensity_range, p0=[halfintensity,x2theta[peak_list[i]],gauss_stddev,halfintensity,lorentz_width])
    peak_parameters.append(parameters)
    x=np.linspace(base_l[i],base_r[i],10000)
    ax.plot(x, _1Voigt(x, *parameters), color="seagreen")
    plt.plot(parameters[1], _1Voigt(parameters[1], *parameters), ".", color="black")
    peak_position.append([parameters[1],np.sqrt(np.diag(info))[1]])

print("-----")
print(filename,"Analysis")
print("-----")

def abs_1Voigt(x, ampG1, cen, sigmaG1, ampL1, widL1):
    return abs((ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen)**2)/((2*sigmaG1)**2)))) +\
              ((ampL1*widL1**2/((x-cen)**2+widL1**2)) ))

peak_intensity=[]
for j in range(len(peak_list)):
    I_j,I_j_info = quad(abs_1Voigt, base_l[j]-peak_width[i]*0.1, base_r[j]+peak_width[i]*0.1, args=(peak_parameters[j][0],peak_parameters[j][1],peak_parameters[j][2],peak_parameters[j][3],peak_parameters[j][4]))
    peak_intensity.append(I_j)
I_max=np.sort(peak_intensity)[len(peak_intensity)-1]
for i in range(len(peak_intensity)):
    peak_intensity[i]/=I_max

print("Normalised peak intensities:",peak_intensity)
print("-----")

plt.show()
