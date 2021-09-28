import numpy as np
import matplotlib.pyplot as plt
import math
from numpy.core.function_base import linspace
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths
import pathlib

#settings

filename = "cu2o_A-X_Y.txt"
cutoff_2xtheta = 12 #degrees
peak_prominence = 20

wavelength = 0.71073 #A^0
max_laue_sum_of_sq = 100

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
plt.xlabel("2x Theta (degrees)")
plt.ylabel("Intensity (arbitrary units)")
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

print("Peak positions and std. dev. (degrees):",peak_position) #entry n corresponds to peak n (n in N0), contains x2theta and std.dev. of peak, respectively. 

sum_of_sq=[]
for i in range(0,max_laue_sum_of_sq+1):
    for j in range(0,max_laue_sum_of_sq+1):
        for k in range(0,max_laue_sum_of_sq+1):
            if i**2+j**2+k**2 not in sum_of_sq and i**2+j**2+k**2<max_laue_sum_of_sq and i**2+j**2+k**2!=0:
                sum_of_sq.append(i**2+j**2+k**2)
sum_of_sq.sort()

cell_parameter_possible=[]
for i in range(len(peak_position)):
    cell_parameter_possible.append([])
    for j in range(len(sum_of_sq)):
        cell_parameter_possible[i].append(math.sqrt(sum_of_sq[j]*wavelength**2 /(4*math.sin(peak_position[i][0]*math.pi/180)**2)))

cell_parameter_possible.sort(reverse=True)
for i in range(len(peak_position)):
    for j in range(len(sum_of_sq)):
        if cell_parameter_possible[i][j]/cell_parameter_possible[0][0]<0.9:
            cell_parameter_possible[i][j]=0
        if cell_parameter_possible[i][j]/cell_parameter_possible[len(peak_position)-1][len(sum_of_sq)-1]>1.1:
            cell_parameter_possible[i][j]=0

cell_parameter_neat=cell_parameter_possible
for i in range(len(peak_position)):
    cell_parameter_neat[i] = [j for j in cell_parameter_possible[i] if j != 0]

print("Possible cell lengths (A^0):",cell_parameter_neat)

plt.show()
