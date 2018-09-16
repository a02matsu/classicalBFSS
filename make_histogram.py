#!/usr/bin/env python3
# coding: UTF-8

from scipy import stats  
import numpy as np  
import matplotlib.pylab as plt
from scipy import optimize

nbin=15

eigenvals = np.loadtxt("EIGENS/XM_Eigen_0002",usecols=(-9,-10))
dist = eigenvals[:,1] - eigenvals[:,0] 

def fit_func(x, b):
    return 2 * b * x * np.exp(-b * x**2)

parameter_initial = np.array([10./dist[0]]) # initial parameters
bin_heights, bin_borders, _ = plt.hist(dist, bins=nbin, label='$N=4$, $t=0$, $\Delta=10$, Sample=1000', normed=True)
bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
popt, pcov = optimize.curve_fit(fit_func, bin_centers, bin_heights , p0=parameter_initial)

print( *popt )


x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
plt.plot(x_interval_for_fit, fit_func(x_interval_for_fit, *popt), label='GOE')
plt.legend()

plt.show()

