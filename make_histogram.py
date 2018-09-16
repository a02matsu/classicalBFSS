#!/usr/bin/env python3
# coding: UTF-8

from scipy import stats  
import numpy as np  
import matplotlib.pylab as plt
from scipy import optimize

nbin=30

eigenvals = np.loadtxt("EIGENS/FM_Eigen_0002",usecols=(9,10))
dist = eigenvals[:,1] - eigenvals[:,0] 

def fit_func_GOE(x, b):
    return 2 * b * x * np.exp(-b * x**2)

def fit_func_GUE(x, c):
    return 4 * c**(3/2) / np.sqrt(np.pi) * x**2 * np.exp(-c * x**2)

#def fit_func_GSE(x, c):
#    return 8 * c**(5/2) / (3*np.sqrt(np.pi)) * x**4 * np.exp(-c * x**2)

## Histram
parameter_initial = np.array([10./dist[0]]) # initial parameters
bin_heights, bin_borders, _ = plt.hist(dist, bins=nbin, label='$N=4$, $t=0$, $\Delta=10$, Sample=100', normed=True)
bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

## fit
popt, pcov = optimize.curve_fit(fit_func_GOE, bin_centers, bin_heights , p0=parameter_initial)
popt2, pcov2 = optimize.curve_fit(fit_func_GUE, bin_centers, bin_heights , p0=parameter_initial)
#popt3, pcov3 = optimize.curve_fit(fit_func_GSE, bin_centers, bin_heights , p0=parameter_initial)



#x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
x_interval_for_fit = np.linspace(0.0, bin_borders[-1], 10000)
plt.plot(x_interval_for_fit, fit_func_GOE(x_interval_for_fit, *popt), label='GOE')
plt.plot(x_interval_for_fit, fit_func_GUE(x_interval_for_fit, *popt2), label='GUE')
#plt.plot(x_interval_for_fit, fit_func_GSE(x_interval_for_fit, *popt3), label='GSE')
plt.legend()

plt.show()

