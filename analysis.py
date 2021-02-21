#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:22:05 2021
@author: leonardobossi1
"""

import numpy as np
import transfermatrix as tmm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

# Importing different materials
Ta2O5 = np.loadtxt("Ta2O5.csv", skiprows=1, unpack=True, delimiter=",")
Ta2O5[0] = Ta2O5*10**3
MgF2 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
BK7 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
Au = np.loadtxt('RefractiveIndexINFO.csv', skiprows = 1, delimiter = ',', unpack = True)

# Converting from micrometers to nanometers
Au[0] = Au[0] * 1000

wavelength, n, k= np.loadtxt("BK7.txt", skiprows = 1, unpack = True)
wavelength = np.array(wavelength) # to avoid problems with curvefit make everything into arrays

n = np.array(n)
k = np.array(k)

# Converted C coefficients to nanometers squared.
S_coefficients = list([1.03961212, 0.231792344, 1.01046945, 6.00069867e3, 2.00179144e4 , 1.03560653e8])

#popt, pcov = curve_fit(sellmeier, wavelength, n, p0 = [1,0,1,0,0,100]) # initial guess taken from internet

# Trial wavelength range
trial = np.arange(330, 2500, 0.2)
interplovals = list(map(lambda x: tmm.linterpol(x, wavelength, n), trial))
sellmeiervals = list(map(lambda x: tmm.sellmeier(x, *S_coefficients), trial))

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

        n1 = tmm.complx_n(incominglam, *MgF2)
        n2 = tmm.complx_n(incominglam, *BK7)

        ns = [n1, n2]
        ds = [j]
        
        #Obtaining the reflection and transmission for a range of values
        
        r, t = tmm.TMM(incominglam, i, incomingpol, ns, ds)
        output.append((i,j,r))
        tot_amp.append(r + t)
    
    # Calculating the theoretical optimal parameters    
    d2 = incominglam/(np.real(n1) * 4 * np.cos(i)) # the analytical formula for d given the phase
    ds = [d2]

    r2, t2 = tmm.TMM(incominglam, i, incomingpol, ns, ds )
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
        
        n1 = tmm.complx_n(lam, *Au) # refractive index of gold
        n2 = tmm.complx_n(lam, *BK7) # refractive index of BK7 glass
        print(type(n1))
        
        ns = [n1, n2]
        ds= [d_val]
        
        r, t = tmm.TMM(lam, ang, incomingpol, ns, ds)

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
    
    n1 = tmm.complx_n(fixed_wavelength, *Au) # refractive index of gold
    n2 = tmm.complx_n(fixed_wavelength, *BK7) # refractive index of BK7 glass
        
    ns = [n1, n2]
    ds= [d_val]
    r, t = tmm.TMM(fixed_wavelength, ang, incomingpol, ns, ds)

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

N, r_current, plot = tmm.find_N(0.9999, fixed_wavelength,  50, 50, incangle, polarisation, Ta2O5, MgF2)

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

# Finding the reflectivity of gold at wavelength = 633 nm (task 12)
n_gold = tmm.complx_n(633, *Au)
n_substrate = tmm.complx_n(633, *BK7)
d_opt = [633 / (np.real(n_gold) * 4)  ]

r_gold = tmm.TMM(633, 0, 's', [n_gold, n_substrate], d_opt)

print('reflectivity of gold at 633nm wavelength = '+str(r_gold[0]))

# Investigating different incoming angles
# Creating a stack
N_stack = 2 # Period of layers
d1 = 50
d2 = 50 # thicknesses of the two material layers

n_stack, d_stack = tmm.stacklayers(N_stack, d1, d2, Ta2O5, MgF2)

ang_range = np.arange(0, np.pi/2, 0.01)
r_output = []
t_output = []
angles= [] 

for ang in ang_range:
    r, t = tmm.TMM(fixed_wavelength, ang, polarisation, n_stack, d_stack)
    
    r_output.append(r)
    t_output.append(t)
    angles.append(ang)
    
#Plotting transmission and reflection
plt.figure()
plt.scatter(angles, r_output, label='R coefficient')
plt.xlabel("Angle of incidence (rad)")
plt.ylabel("Reflectance/Transmission coefficient")
plt.legend()
plt.grid()

plt.figure()
plt.scatter(angles, t_output, label='T coefficient', color='r')
plt.xlabel("Angle of incidence (rad)")
plt.ylabel("Reflectance/Transmission coefficient")
plt.legend()
plt.grid()    
    
#%%
# for wavelength of 633nm, phase response of reflected wave with both gold and dbr
n_gold = tmm.complx_n(633, *Au)
n_substrate = tmm.complx_n(633, *BK7)
dG = [633 / (np.real(n_gold) * 4)]

nT = tmm.complx_n(fixed_wavelength,*Ta2O5)
nM = tmm.complx_n(fixed_wavelength,*MgF2)
dT = 633/(np.real(nT)*4*np.cos(incangle))
dM = 633/(np.real(nM)*4*np.cos(incangle))

n_stack, d_stack = tmm.stacklayers(14, 633, dM, dT, MgF2, Ta2O5, nfinal = n_substrate)

var_wavelength = np.arange(551,741)
rsphases = []
rgphases = []
for i in var_wavelength:
    rstack, tstack = tmm.TMM(i, 0, "s", n_stack,d_stack, squared = False)
    rgold, tgold = tmm.TMM(i, 0, "s", [n_gold,n_substrate],dG, squared = False)
    rsphase = np.angle(rstack)
    rgphase = np.angle(rgold)
    rsphases.append(rsphase)
    rgphases.append(rgphase)

plt.figure()
plt.plot(var_wavelength,rsphases, label = "Phase response of DBR")
plt.plot(var_wavelength,rgphases, label = "Phase response of gold")
plt.legend()
plt.xlabel("Wavength in nm")
plt.ylabel("Phase of reflectance in rad")
plt.yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
          [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$'])
plt.show()   
    





    
    
    
   
    




