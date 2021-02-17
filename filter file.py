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

#%%
def chifunction (k0, kx, polarisation, ncomplex1, ncomplex2):

    n1 = np.real(ncomplex1)
    n2 = np.real(ncomplex2)

    kz1 = np.sqrt((n1*k0)**2 - kx**2)
    kz2 = np.sqrt((n2*k0)**2 - kx**2)

    if polarisation == "p":

        alpha = (n2/n1)*np.sqrt(kz1/kz2)
        chip = (alpha + 1/alpha)/2
        chim = (alpha - 1/alpha)/2

    elif polarisation == "s":

        alpha = np.sqrt(kz2/kz1)
        chip = (alpha + 1/alpha)/2
        chim = (alpha - 1/alpha)/2

    else:
        raise Exception("That is not a supported polarisation of light")
    
    return chip, chim # returns chi plus and chi minuts    
    
def TMM(wavelength, angle, polarisation, ns, ds):

    if polarisation != 's' and polarisation != 'p':
        raise Exception("That is not a supported polarisation of light")

    k0 = 2*np.pi/wavelength
    kx = k0*np.sin(angle)
    M = [[1,0],[0,1]] # initialise general matrix
    for i in range(len(ds)): # ds should be one item shorter than ns, the final ns should be for the substrate
        n = np.real(ns[i])
        k = np.imag(ns[i])
        kz = np.sqrt((n*k0)**2 - kx**2)

        if i == 0:
            n1 = 1 # air
            n2 = ns[i] # the refractive index of the layer

        else:
            n1 = ns[i-1]
            n2 = ns[i]

        propagation = np.exp(complex(0,(kz*ds[i]))) # forward propagation
        chip, chim = chifunction(k0, kx, polarisation, n1, n2)
        T_i = [[chip , chim],[chim, chip]]
        P_i = [[propagation, 0],[0, np.conj(propagation)]] # complex conjugate for backward propagation
        interstep = np.matmul(P_i, T_i)
        M = np.matmul(interstep, M)

    n1 = ns[-2]
    n2 = ns[-1]

    chip, chim = chifunction(k0, kx, polarisation, n1, n2)
    T_i = [[chip , chim],[chim, chip]] # interfacial for the substrate
    M = np.matmul(T_i, M)

    r = -(M[1][0]/M[1][1])
    t = M[0][0] + M[0][1]*r
    
    r = r**2
    t = t**2 # want the transmittance and reflectance
    
    return r, t
  
data2 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
data1 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
incominglam = 500
d = np.arange(0, 160, 0.25)
incangle = np.arange(0, 1, 0.01)
incomingpol = "s"

output = []
analytic = []



for i in incangle:
    for j in d: # numerically getting all values

        n1 = complx_n(incominglam, *data1)
        n2 = complx_n(incominglam, *data2)

        ns = [n1, n2]
        ds = [j]

        r, t = TMM(incominglam, i, incomingpol, ns, ds)
        output.append((i,j,r))
    d2 = incominglam/(np.real(n1)*4*np.cos(i)) # the analytical formula for d given the phase
    ds = [d2]

    r2, t2 = TMM(incominglam, i, incomingpol, ns, ds )
    analytic.append((i, d2, r2))

xcoord = []
ycoord = []
zcoord = []
for i in output:

    xcoord.append(i[0])
    ycoord.append(i[1])
    zcoord.append(i[2])

xcoord2, ycoord2, zcoord2 = [], [], []
for i in analytic:

    xcoord2.append(i[0])
    ycoord2.append(i[1])
    zcoord2.append(i[2])

ycoords = ycoord[::4]
xcoords = xcoord[::4]
zcoords = zcoord[::4]

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel("incidence angle rad")
ax.set_ylabel("thickness nm")
ax.set_zlabel("reflection")


ax.plot(xcoord2, ycoord2, zcoord2, color = 'r', zorder = 1, label =  "analytical formula")
ax.scatter(xcoord, ycoord, zcoord, c=zcoord,cmap=cm.viridis, zorder = 2, alpha = 0.05)

ax.legend()

plt.show()

















