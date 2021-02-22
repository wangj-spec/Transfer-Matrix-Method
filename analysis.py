#!/usr/bin/env python3
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
MgF2 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
BK7 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
Au = np.loadtxt('RefractiveIndexINFO.csv', skiprows = 1, delimiter = ',', unpack = True)

Ta2O5[0] = Ta2O5[0] * 1000

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
visible_spec = np.arange(500, 900, 1)
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
import transfermatrix as tmm

fixed_wavelength = 633
incangle = 0
polarisation = "s"

# Finding the number of periods required to reach 99.99% reflectivity, d1=d2=50nm,
# We used materials Ta2O5 and MgF2 with air substrate.

lam_opt = 633

d1opt =lam_opt / (4 * np.real(tmm.complx_n(lam_opt, *Ta2O5)))
d2opt = lam_opt / (4 * np.real(tmm.complx_n(lam_opt, *MgF2)))

N, r_current, plot = tmm.find_N(0.9999, fixed_wavelength,  d1opt, d2opt, incangle, polarisation, Ta2O5, MgF2)
nplot = []
rplot = []

for i in plot:
    nplot.append(i[0])
    rplot.append(i[1])

im = plt.figure()
plt.xlabel("Number of stacks")
plt.ylabel("Reflectance")
plt.title('Number of stacks (each stack has 2 layers) required to reach 99.99% reflectivity.')
plt.scatter(nplot, rplot)
plt.grid()
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

n_stack, d_stack = tmm.stacklayers(N_stack, 500, d1opt, d2opt, Ta2O5, MgF2)

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


# Trying to find the stop band by varying wavelength
# We use the optimal thicknesses for an incoming wavelength of 633nm
# so this is where the peak should be observed.

N_range = np.arange(1, 10, 2)
incangle = 0

for N in N_range:
    r_values = []    
    
    for lam in visible_spec:
        n_stack2, d_stack2 = tmm.stacklayers(N, lam, d1opt, d2opt, Ta2O5, MgF2)
        
        r, t = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2)
        r_values.append(r)
    
    
    plt.plot(visible_spec, r_values, label='N = '+str(N))

plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient ')
plt.title('Reflectivity spectrum as a function of wavelength for different N values')
plt.legend(loc='best')
plt.grid()
        
    
# Plotting spectrum from N=14 (layers required to reach 99.99% reflectivity).
# Can also be used to investigate band-width change for different materials by
# changing "material 1 and material 2)
plt.figure()

N = 12
r_values = []

layer1 = Ta2O5
layer2 = MgF2

# Ensuring the thicknesses are optimal for reflectivity 
d1 =lam_opt / (4 * np.real(tmm.complx_n(lam_opt, *layer1)))
d2 = lam_opt / (4 * np.real(tmm.complx_n(lam_opt, *layer2)))

for lam in visible_spec:
    n_stack2, d_stack2 = tmm.stacklayers(N, lam, d1, d2, layer1, layer2)
        
    r, t = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2)
    r_values.append(r)
    
plt.plot(visible_spec, r_values, label='N = '+str(N))
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient')
plt.legend(loc='best')
plt.grid()
        

import matplotlib as mpl
# Set the default color cycle
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["r", "k", "c"]) 
    
angles = [0 ,np.pi/6, np.pi/3]

for ang in angles:
    
    r_values = []
    
    for lam in visible_spec:
        n_stack2, d_stack2 = tmm.stacklayers(14, lam, d1opt, d2opt, Ta2O5, MgF2)
        
        r, t = tmm.TMM(lam, ang, polarisation, n_stack2, d_stack2)
        r_values.append(r)
    
    if ang !=0:        
        plt.plot(visible_spec, r_values, label='Incident Angle = \u03C0 /'+str(int(1/ang*(np.pi))))
    else:
        plt.plot(visible_spec, r_values, label='Incident Angle = 0', color='b', linestyle='--')
    
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient ')
plt.title('Reflectivity spectrum as a function of wavelength for different incident angles')
plt.legend(loc='best')
plt.grid()
    

# Obtaining the plot for the extreme case, (as incident angle tends towards pi/2).
# Stacks used 14
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["black", "dimgray", "r"]) 
angle_tolim = [2 * np.pi/ 6, 11 * np.pi / 12, 99 * np.pi/ 200]

for ang in angle_tolim:
    r_values = []
    for lam in visible_spec:
        
        n_stack2, d_stack2 = tmm.stacklayers(14, lam, d1opt, d2opt, Ta2O5, MgF2)
            
        r, t = tmm.TMM(lam, ang, polarisation, n_stack2, d_stack2)
        r_values.append(r)
        
    plt.plot(visible_spec, r_values)
    
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient ')
plt.title('Reflectivity as angle tends towards \u03C0 / 2 ')
plt.legend(loc='best')
plt.grid()

# Investigating changing materials and the effect it has on spectral width

def spectral_width(r_data):
    '''
    Parameters:
        r_data::list
            Reflection coefficient. (for increasing wavelengths.)
    Returns:
        Width::float
            Spectral width at FWHM
    '''
    
    
    lower_index = r_data.index(max(r_data))
    upper_index = r_data.index(max(r_data))

    while r_data[upper_index] > 0.5:
        upper_index+=1
    while r_data[lower_index] > 0.5:
        lower_index-=1
    
    width = visible_spec[upper_index] - visible_spec[lower_index]
    
    return width
    
    
# Using a stack with different materials (SiO2 with MgF2, TiO2 with Ta2O5)
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["r", "k", "c"]) 

TiO2 = np.loadtxt("Devore-o.csv", skiprows=1, unpack=True, delimiter=",")
SiO2 = np.loadtxt("Malitson.csv", skiprows=1, unpack=True, delimiter=",")

TiO2[0] = TiO2[0] * 1000
SiO2[0] = SiO2[0] * 1000

dopt_TiO2 = lam_opt / (4 * np.real(tmm.complx_n(lam_opt, *TiO2)))
dopt_SiO2 = lam_opt / (4 * np.real(tmm.complx_n(lam_opt, *SiO2)))

N = 14

# TiO2 combined with Ta2O5 (both high refractive indices) and TiO2 with SiO2 
# (very similar refractive indices as materials we already used)
plt.figure()

r_values = []
r_values2 = []

for lam in visible_spec:
    n_stack, d_stack = tmm.stacklayers(N, lam, d1opt, dopt_TiO2, Ta2O5, TiO2)
    r, t = tmm.TMM(lam, incangle, polarisation, n_stack, d_stack)
    r_values.append(r)
    
    n_stack2, d_stack2 = tmm.stacklayers(N, lam, dopt_TiO2, dopt_SiO2, TiO2, SiO2)
    r2, t2 = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2)
    r_values2.append(r2)
    

plt.plot(visible_spec, r_values, label='N=14 using TiO2 and Ta2O5')
plt.plot(visible_spec, r_values2, label='N =14 using TiO2 and SiO2')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient ')
plt.grid()
plt.legend(loc='upper right')

# MgF2 and SiO2 (both lower refractive indices). Many more layers are required
# to achieve high reflectivity, and the spectral width is much lower.

N= 40
r_values3 = []

for lam in visible_spec:
    n_stack2, d_stack2 = tmm.stacklayers(N, lam, dopt_SiO2, d2opt, SiO2, MgF2)
        
    r, t = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2)
    r_values3.append(r)

plt.figure()
plt.plot(visible_spec, r_values3, label='N=40 using SiO2 and MgF2', color='r')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient')
plt.grid()
plt.legend()
    
print('spectral width for SiO2 and MgF2 with N=35 is '+str(spectral_width(r_values3))+' nm')
print('spectral width for TiO2 with SiO2 with N=14 is '+str(spectral_width(r_values2))+' nm')
print('spectral width for TiO2 and Ta2O5 with N=14 is '+str(spectral_width(r_values))+' nm')
