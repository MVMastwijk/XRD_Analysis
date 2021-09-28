# XRD Powder Diffraction Pattern Analysis Program (Cubic Crystal Sytems)

the data file must be in the same folder as the program


settings:

filename: name of data file

cutoff_x2theta: (2x) angle from which peaks are counted

peak_prominence: required prominence of peaks (high value means only tallest peaks get selected)

wavelength: wavelength of x-rays

laue_indices_sum_of_sq: list of sum of squares of laue indices used in analysis of the data


output:

plot of xrd data including gaussian fit and maxima of fit

all (2x) diffraction angles + standard deviations for the peaks (calculated by single gaussian function fitting over each individual peak)

the values of 4sin(θ)^2/λ^2, used to determine the laue indices sum of squares 

the cell parameter and volume, including standard deviations
