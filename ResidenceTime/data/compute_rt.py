#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp

def homogenous_poisson_process_cdf(t, rate):
    # rate = 1/tau
    return 1-np.exp(-t*rate)

data = np.loadtxt("./simulations_data.dat")

range_factor = 100
nbins = 10_000
bins_time = np.logspace(
    np.log10(np.amin(data)/range_factor),
    np.log10(np.amax(data)*range_factor),
    num=nbins, base=10
    )
hist_values, _ = np.histogram(data, bins=bins_time)
hist_values = np.append(hist_values,0)
cdf_e = np.cumsum(hist_values)/len(data)

popt, _ = curve_fit(
    homogenous_poisson_process_cdf,
    xdata= bins_time,
    ydata=cdf_e,
    p0=1/np.mean(data))
tau = 1/popt[0]

sampling_from_theoretical_distribution = np.random.exponential(
    scale=tau,
    size=int(len(data)*1E6)
    )
_, ks_p_value = ks_2samp(data, sampling_from_theoretical_distribution)


print(f"Residence time from fitting: {tau:.4e} [in the units of the original data]")
print(f"p-value: {ks_p_value:.3f}")
