#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:22:05 2021

@author: leonardobossi1
"""

import scipy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linterpol(wavelength, wavedata, indexdata):

    grad = 0 #initialise grad
    j = 0 #initialise j

    for i in range(len(wavedata)):
        if wavedata[0] >= wavelength:
            grad = (wavelength - wavedata[0])/(wavedata[1] - wavedata[0])
            j = 1
            break # linear extrapolation
        elif wavedata[i] >= wavelength:
            grad = (wavelength - wavedata[i])/(wavedata[i] - wavedata[i-1])
            j = i
            break # linear interpolation
        elif wavedata[-1] < wavelength:
            grad = (wavelength - wavedata[-1])/(wavedata[-1] - wavedata[-2])
            j = -1
            break # linear extrapolation

    interpol = grad*(indexdata[j] - indexdata[j-1]) + indexdata[j]
    
    return interpol


def sellmeier(wavelength, a1, a2, a3, b1, b2, b3):
    n_squared = 1 + (a1 * wavelength ** 2)/(wavelength ** 2 - b1) + (a2 * wavelength **2 )/ (wavelength ** 2 - b2) +\
             (a3 * wavelength ** 2)/(wavelength ** 2 - b3)

    n = np.sqrt(n_squared)
    
    return n


wavelength, n, k= np.loadtxt("BK7.txt", skiprows = 1, unpack = True)
wavelength = np.array(wavelength) # to avoid problems with curvefit make everything into arrays

n = np.array(n)
k = np.array(k)

# Converted C coefficients to nanometers squared.
S_coefficients = list([1.03961212, 0.231792344, 1.01046945, 6.00069867e3, 2.00179144e4 , 1.03560653e8])

#popt, pcov = curve_fit(sellmeier, wavelength, n, p0 = [1,0,1,0,0,100]) # initial guess taken from internet

# Trial wavelength range
trial = np.arange(330, 2500, 0.2)
interplovals = list(map(lambda x: linterpol(x, wavelength, n), trial))
sellmeiervals = list(map(lambda x: sellmeier(x, *S_coefficients), trial))

plt.figure()
plt.plot(trial, interplovals)
plt.plot(trial, sellmeiervals)
plt.grid()
residuals = []

for i in range(len(trial)):
    residuals.append(interplovals[i] - sellmeiervals[i])

print("The averaage residual is " + str(np.mean(residuals)))

plt.figure()
plt.scatter(trial, residuals)
plt.title('Residual values between using linear interpolating and Sellmeier equation')
plt.grid()

plt.figure()
plt.title('Scatter plot of real values from refractive indices')
plt.plot(wavelength, n, 'x', label = 'real values for BK7 glass')
plt.xlabel('Wavelength (nm)')
plt.ylabel('refractive index')
plt.grid()
plt.legend()

plt.figure()
plt.title('Scatter plot of imaginary values from refractive indices')
plt.xlabel('Wavelength (nm)')
plt.ylabel('refractive index')
plt.plot(wavelength, k, 'x', label = 'imaginary values for BK7 glass')
plt.grid()

#%%
# Gold data, wavelength is given in micrometers

gold_data = np.loadtxt('RefractiveIndexINFO.csv', skiprows = 1, delimiter = ',', unpack = True)

# Converting to nonmeters
Au_wavelgth = gold_data[0] * 10 ** 3 


def complx_n(lam, lam_data = Au_wavelgth, real_data = gold_data[1], img_data = gold_data[2]):
    
    if lam > lam_data[-1] or lam < lam_data[0]:
        raise Exception('Inputted value is out of range provided in the data')
        
    n_real = linterpol(lam, lam_data, real_data)
    n_img = linterpol(lam, lam_data, img_data)
    n = complex(n_real, n_img)
    
    return n

    
    
    
















