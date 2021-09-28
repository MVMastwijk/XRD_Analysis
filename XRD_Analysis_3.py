#XRD Powder Diffraction Analysis Program (for cubic crystal systems)
#Author = MV Mastwijk, 28-SEP-2021
#See below for the settings

import numpy as np
import matplotlib.pyplot as plt
import math
from numpy.core.function_base import linspace
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths
import pathlib

#settings

filename = "cu2o_A-X_Y.txt" #note: file must be in the same folder as the program
cutoff_2xtheta = 12 #degrees
peak_prominence = 20 #higher peak prominence > small peaks are ignored

wavelength = 0.71073 #A^0
laue_indices_sum_of_sq = [2,3,4,6,8,11,12,16,19,20] #when first running program, leave empty. determine the values by hand via 4sin(θ)^2/λ^2=l/a^2 and rerun

#

powder_diffraction = np.genfromtxt(pathlib.PurePath(pathlib.Path(__file__).parent, filename))

x2theta,intensity=[],[]
for i in range(len(powder_diffraction)):
    x2theta.append(powder_diffraction[i][0])
    intensity.append(powder_diffraction[i][1])

fig, ax = plt.subplots(1, figsize=(8, 6), facecolor='whitesmoke')
ax.set_facecolor("gainsboro")
ax.grid(color='white', linestyle='-', linewidth=1)
fig.suptitle(filename, fontsize=15)
plt.xlabel("2θ (°)")
plt.ylabel("I (arbitrary units)")
ax.plot(x2theta, intensity, color="orange")
plt.axvline(x=cutoff_2xtheta, linestyle='--', color="grey")

peak_list, peak_data=find_peaks(intensity, prominence=peak_prominence, width=0, rel_height=1)
peak_base_height=peak_data.get("width_heights")
peak_width=peak_widths(intensity,peak_list,rel_height=0.5)[0]

peak_x,peak_y=[],[]
base_l,base_r=[],[]
n=0
for i in range(len(peak_list)):
    if x2theta[peak_list[i]]<cutoff_2xtheta:
        n+=1
    else:
        peak_x.append(x2theta[peak_list[i]])
        peak_y.append(intensity[peak_list[i]])
        peak_width[i]=int(round(peak_width[i]))
        base_l.append(x2theta[int(peak_list[i]-peak_width[i])])
        base_r.append(x2theta[int(peak_list[i]+peak_width[i])])
for i in range(n):
    peak_list=np.delete(peak_list, 0)
    peak_base_height=np.delete(peak_base_height, 0)
    peak_width=np.delete(peak_width, 0)

def gaussian(x, amplitude, mean, stddev):
    return baseheight+amplitude * np.exp(-((x - mean) / 4 / stddev)**2)

peak_position=[]
for i in range(len(peak_list)):
    intensity_range=[]
    x2theta_range=[]
    baseheight=peak_base_height[i]
    for j in range(len(x2theta)):
        if x2theta[j]>=base_l[i] and x2theta[j]<=base_r[i]:
            x2theta_range.append(x2theta[j])
            intensity_range.append(intensity[j])
    parameters, info = curve_fit(gaussian, x2theta_range, intensity_range, p0=[intensity[peak_list[i]],x2theta[peak_list[i]],1])
    x=np.linspace(base_l[i],base_r[i],10000)
    ax.plot(x, gaussian(x, *parameters), color="seagreen")
    plt.plot(parameters[1], parameters[0]+baseheight, ".", color="black")
    peak_position.append([parameters[1],np.sqrt(np.diag(info))[1]])

print("-----")
print("Peak positions and std. dev. (degrees):",peak_position) #entry n corresponds to peak n (n in N0), contains x2theta and std.dev. of peak, respectively. 
print("-----")

factor_list=[]
for i in range(len(peak_position)):
    factor_list.append((4*math.sin(0.5*peak_position[i][0]*math.pi/180)**2)/wavelength**2)

cell_parameter=[]
cell_parameter_variance=[]
cell_volume=[]
cell_volume_variance=[]
firstrun=False
if laue_indices_sum_of_sq==[]:
    for i in range(len(peak_position)):
        laue_indices_sum_of_sq.append(1)
        firstrun=True

for i in range(len(peak_position)):
    cell_ai=1/math.sqrt(factor_list[i]/laue_indices_sum_of_sq[i])
    cell_parameter.append(cell_ai)
    cell_parameter_variance.append((cell_ai**2)/(math.tan(0.5*peak_position[i][0]*math.pi/180)**2)*(0.5*peak_position[i][1]*math.pi/180)**2)
    cell_volume.append(cell_ai**3)
    cell_volume_variance.append(9*cell_ai**4*cell_parameter_variance[i])

cell_parameter_avg=sum(cell_parameter)/len(peak_position)
cell_parameter_variance_avg=sum(cell_parameter_variance)/(len(peak_position)**2)

cell_volume_avg=sum(cell_volume)/len(peak_position)
cell_volume_variance_avg=sum(cell_volume_variance)/(len(peak_position)**2)


print("4sin(θ)^2/λ^2:",factor_list)
print("-----")
if firstrun==False:
    print("Sum of sq. of Laue indices:",laue_indices_sum_of_sq)
    print("-----")
    print("Cell parameter (A^0):",cell_parameter_avg)
    print("Cell parameter std. dev. (A^0):",math.sqrt(cell_parameter_variance_avg))
    print("-----")
    print("Cell volume (A^0^3):",cell_volume_avg)
    print("Cell volume std. dev. (A^0^3):",math.sqrt(cell_volume_variance_avg))
else:
    print("Determine the sum of squares l of Laue indices belonging to each peak by 4sin(θ)^2/λ^2=l/a^2, put these in the settings in ascending order and rerun.")
print("-----")

plt.show()
