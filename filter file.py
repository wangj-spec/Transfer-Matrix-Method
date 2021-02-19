#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:22:05 2021
@author: leonardobossi1
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

# Importing different materials
Ta2O5 = np.loadtxt("Ta2O5.csv", skiprows=1, unpack=True, delimiter=",")
MgF2 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
BK7 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
Au = np.loadtxt('RefractiveIndexINFO.csv', skiprows = 1, delimiter = ',', unpack = True)

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


def complx_n(lam, lam_data, real_data, img_data):
    
    if lam > lam_data[-1] or lam < lam_data[0]:
        raise Exception('Inputted value is out of range provided in the data')
        
    n_real = linterpol(lam, lam_data, real_data)
    n_img = linterpol(lam, lam_data, img_data)
    n = complex(n_real, n_img)
    
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

print("The average residual is " + str(np.mean(residuals)))

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

    k0 = 2 * np.pi / wavelength
    kx = k0 * np.sin(angle)
    M = [[1,0],[0,1]] # initialise general matrix
    
    for i in range(len(ds)): # ds should be one item shorter than ns, the final ns should be for the substrate
        n = np.real(ns[i])
        k = np.imag(ns[i])  
        kz = np.sqrt((n * k0) ** 2 - kx ** 2)
        absorp_factor = - k * ds[i] * kz / n # Absorption coefficient in layer 

        if i == 0: # setting incident medium to be air
            n1 = 1 # air
            n2 = ns[i] # the refractive index of the layer

        else:
            n1 = ns[i-1]
            n2 = ns[i]

        propagation = np.exp(complex(absorp_factor,(kz * ds[i]))) # forward propagation
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
    
    r = abs(r) ** 2
    t = abs(t) ** 2 # want the transmittance and reflectance
    
    return r, t
  
    
# incominglam is the wavelength of the incident light
# d is the thickness of the anti-reflection layer
# incangle is the angle of incidence

incominglam = 500
d = np.arange(0, 160, 0.25)
incangle = np.arange(0, 1, 0.01)
incomingpol = "s"

output = []
analytic = []
tot_amp = []

for i in incangle:
    for j in d: # numerically getting all values

        n1 = complx_n(incominglam, *MgF2)
        n2 = complx_n(incominglam, *BK7)

        ns = [n1, n2]
        ds = [j]
        
        #Obtaining the reflection and transmission for a range of values
        
        r, t = TMM(incominglam, i, incomingpol, ns, ds)
        output.append((i,j,r))
        tot_amp.append(r + t)
    
    # Calculating the theoretical optimal parameters    
    d2 = incominglam/(np.real(n1) * 4 * np.cos(i)) # the analytical formula for d given the phase
    ds = [d2]

    r2, t2 = TMM(incominglam, i, incomingpol, ns, ds )
    analytic.append((i, d2, r2))

xcoord = []
ycoord = []
zcoord = []

for i in output: # Output consists of angle value, thickness value, and reflection coefficient

    xcoord.append(i[0])
    ycoord.append(i[1])
    zcoord.append(i[2])

# Creating the array of analytical values expected
    
xcoord2, ycoord2, zcoord2 = [], [], []

for i in analytic:

    xcoord2.append(i[0])
    ycoord2.append(i[1])
    zcoord2.append(i[2])

ycoords = ycoord[::10]
xcoords = xcoord[::10]
zcoords = zcoord[::10]

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel("incidence angle rad")
ax.set_ylabel("thickness nm")
ax.set_zlabel("reflection")


ax.plot(xcoord2, ycoord2, zcoord2, color = 'r', zorder = 1, label =  "analytical formula")
ax.scatter(xcoord, ycoord, zcoord, c=zcoord, cmap=cm.viridis, zorder = 2, alpha = 0.05)

ax.legend()

plt.show()

#  Checking if r + t = 1 (it should be equal to 1 as we haven't considered
# absorption yet).

print('Maximum t+r value = '+str(max(tot_amp))+', minimum value ='+str(min(tot_amp))+ \
      ' Both of these values are on the order of 10^-16, suggesting this this a floating '\
      'point error.')



 
#%%
# Absorption for a layer of gold
    
# Investigating for normal incidience (angle = 0)

ang = 0 
visible_spec = np.arange(380, 750, 1)
d_range = np.arange(10, 100, 1)
incomingpol = "s"

# Creating output list with: wavelength, thickness, reflection, and transmission
output_list = []

for lam in visible_spec:
    for d_val in d_range:
        
        n1 = complx_n(lam, *Au) # refractive index of gold
        n2 = complx_n(lam, *BK7) # refractive index of BK7 glass
        
        ns = [n1, n2]
        ds= [d_val]
        
        r, t = TMM(lam, ang, incomingpol, ns, ds)

        output_list.append((lam, d_val, r, t))
        
    
xcoord2, ycoord2, r_vals, t_vals = [], [], [], []

for output in output_list:

    xcoord2.append(output[0])
    ycoord2.append(output[1])
    r_vals.append(output[2])
    t_vals.append(output[3])


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel("Wavelegth")
ax.set_ylabel("thickness nm")
ax.set_zlabel("reflection")
ax.set_title("Reflection variation for Gold layer on top of glass substrate for normal incidence")


ax.scatter(xcoord2, ycoord2, r_vals, c = t_vals, cmap=cm.viridis, zorder = 2, alpha = 0.05)

ax.legend()

# Plotting the transmssion spectrum as a function of wavelength and thickness
fig2 = plt.figure()
ax2 = plt.axes(projection='3d')

ax2 = plt.axes(projection='3d')
ax2.set_xlabel("Wavelegth")
ax2.set_ylabel("thickness nm")
ax2.set_zlabel("Transmission coefficient")
ax2.set_title("Transmission spectrum for Gold layer on top of glass substrate for normal incidence")


ax2.scatter(xcoord2, ycoord2, t_vals, c=t_vals, cmap=cm.viridis, zorder = 2, alpha = 0.05)

ax2.legend()


data2d = []


fixed_wavelength = 500
ang = 0 # normal incidence

for d_val in d_range:
    
    n1 = complx_n(fixed_wavelength, *Au) # refractive index of gold
    n2 = complx_n(fixed_wavelength, *BK7) # refractive index of BK7 glass
        
    ns = [n1, n2]
    ds= [d_val]
    r, t = TMM(fixed_wavelength, ang, incomingpol, ns, ds)

    data2d.append((d_val, r, t))
    
data2d = np.array(data2d)
    
plt.figure()
plt.scatter(data2d[:,0], data2d[:,2], label = 'Transmission spectrum for fixed wavelength 500nm')
plt.title('Gold layer with glass substrate')
plt.xlabel('Thickness of layer (nm)')
plt.ylabel('Transmission coefficient')
plt.grid()
        
 
#%%

fixed_wavelength = 633
incangle = 0
polarisation = "s"

def stacklayers(N, d1 , d2, material1, material2):
    ns = []
    ds = []
    n1 = complx_n(fixed_wavelength,*material1)
    n2 = complx_n(fixed_wavelength,*material2)


    for i in range(N):
        ns.append(n1)
        ds.append(d1)
        ns.append(n2)
        ds.append(d2)

    return ns, ds



def find_N(r_val, wavelength, d1, d2, angle, polarisation, material1, material2):
    N = 1
    plot = []
    r_current = 0 # initialies r_current
    
    while r_current < r_val:
        ns, ds = stacklayers(N, d1, d2, material1, material2)
        r_current = TMM(wavelength, angle, polarisation, ns, ds)[0]
        plot.append([N, r_current])
        N += 1
    plot.append([N, r_current])

    return N, r_current, plot

n_1 = complx_n(fixed_wavelength, *MgF2)
n_1 = np.real(n_1)
n_2 = complx_n(fixed_wavelength, *Ta2O5)
n_2 = np.real(n_2)

N, r_current, plot = find_N(0.9999, fixed_wavelength,  50, 50, incangle, polarisation, Ta2O5, MgF2)

nplot = []
rplot = []

for i in plot:
    nplot.append(i[0])
    rplot.append(i[1])

im = plt.figure()
plt.xlabel("Number of stacks")
plt.ylabel("Reflectance")
plt.scatter(nplot, rplot)
plt.show() 
        

# Investigating different incoming angles

# Creating a stack
N_stack = 2 # Period of layers
d1 = 50
d2 = 50 # thicknesses of the two material layers

n_stack, d_stack = stacklayers(N_stack, d1, d2, Ta2O5, MgF2)

ang_range = np.arange(0, np.pi/2, 0.01)
r_output = []
t_output = []
test= [] 

for ang in ang_range:
    r, t = TMM(fixed_wavelength, ang, polarisation, n_stack, d_stack)
    
    r_output.append(r)
    t_output.append(t)
    test.append(ang)
    
    
    
   
    




