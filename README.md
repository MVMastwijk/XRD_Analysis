# XRD_Analysis (not finalized)

data file must be in the same folder as the program
!note the program assumes a cubic crystal system


settings:

filename: name of data file

cutoff_x2theta: (2x) angle from which peaks are counted

peak_prominence: required prominence of peaks (high value means only tallest peaks get selected)

wavelength: wavelength of x-rays

max_laue_sum_of_sq: maximal sum of squares of laue indices used in analysis of the data


output:

plot of xrd data including gaussian fit and maxima of fit

all (2x) diffraction angles + standard deviations for the peaks (calculated by single gaussian function fitting over each individual peak)

possible cell parameters (for chosen max_laue_sum_of_sq)
