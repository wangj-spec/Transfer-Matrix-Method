"""
Created on Wed Mar  3 00:43:27 2021

@author: leonardobossi1
"""

import numpy as np
import tmmfile as tmm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

# Importing different materials
MgF2 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
BK7 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
Au = np.loadtxt('Au.txt', skiprows=1, unpack=True)

wavelength, n, k = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
wavelength = np.array(wavelength)  # to avoid problems with curvefit make everything into arrays

n = np.array(n)
k = np.array(k)

# Converted C coefficients to nanometers squared.
S_coefficients = list([1.03961212, 0.231792344, 1.01046945, 6.00069867e3, 2.00179144e4, 1.03560653e8])


# Trial wavelength range
trial = np.arange(330.5, 2500, 0.2)
interplovals = list(map(lambda x: tmm.linterpol(x, wavelength, n), trial))
sellmeiervals = list(map(lambda x: tmm.sellmeier(x, *S_coefficients), trial))

plt.figure()
plt.title("Refractive index of BK7")
plt.plot(trial, interplovals)
plt.plot(trial, sellmeiervals)
plt.xlabel('Wavelength (nm)')
plt.ylabel("refractive index")
plt.grid()
residuals = []

for i in range(len(trial)):
    residuals.append(interplovals[i] - sellmeiervals[i])

print("The average residual is " + str(np.mean(residuals)))

plt.figure()
plt.scatter(trial, residuals)
plt.title('Residual values between using linear interpolating and Sellmeier equation')
plt.xlabel('Wavelength (nm)')
plt.ylabel("Residual values")
plt.grid()

plt.figure()
plt.title('Scatter plot of real values from refractive indices')
plt.plot(wavelength, n, 'x', label='real values for BK7 glass')
plt.xlabel('Wavelength (nm)')
plt.ylabel('refractive index')
plt.grid()
plt.legend()

plt.figure()
plt.title('Scatter plot of imaginary values from refractive indices')
plt.xlabel('Wavelength (nm)')
plt.ylabel('refractive index')
plt.plot(wavelength, k, 'bx', label='imaginary values for BK7 glass')
plt.grid()

# %%
# incominglam is the wavelength of the incident light
# d is the thickness of the anti-reflection layer
# incangle is the angle of incidence

incominglam = 500 # Fixing the wavelength to 500nm and seeing how the reflectance
                  # varies with the other variables.

d = np.arange(0, 160, 0.25)
incangle = np.arange(0, 1, 0.01)
incomingpol = "p"

plt.figure()
#To see the reflectance varying with distance
r_values = []

test = np.arange(0, 2, 0.01)
n1 = tmm.complx_n(incominglam, *MgF2)
n2 = tmm.complx_n(incominglam, *BK7)
ns = [n1, n2]

for i in test:
    dm = i*np.real(incominglam / (2 * (tmm.complx_n(incominglam, *MgF2))))
    ds = [dm]
    r, t = tmm.TMM(incominglam, 0, incomingpol, ns, ds)
    r_values.append(r)

plt.plot(test, r_values)

plt.xlabel('Phase due to width of layer (*pi)')
plt.ylabel('Reflection coefficient ')
plt.title('Reflectivity spectrum as a function of phase')
plt.legend(loc='best')
plt.grid()
plt.show()

output = []
analytic = []
tot_amp = []

for i in incangle:
    for j in d:  # numerically getting all values

        n1 = tmm.complx_n(incominglam, *MgF2)
        n2 = tmm.complx_n(incominglam, *BK7)
        n1 = complex(np.real(n1), 0)
        n2 = complex(np.real(n2), 0)
        ns = [n1, n2]
        ds = [j]

        # Obtaining the reflection and transmission for a range of values

        r, t, a = tmm.TMM(incominglam, i, incomingpol, ns, ds, absorption=True)
        output.append((i, j, r))
        tot_amp.append(r + t)

    # Calculating the theoretical optimal parameters for each angle
        ko = 2 * np.pi/incominglam
        kx = ko * np.sin(i)
        
        d2 = np.pi/ (2 * np.sqrt((np.real(n1)*ko)**2 - kx**2))  # the analytical formula for d given the phase
        #print(np.shape(d2))
        ds = [d2]

        r2, t2 = tmm.TMM(incominglam, i, incomingpol, ns, ds)
        analytic.append((i, d2, r2))

xcoord = []
ycoord = []
zcoord = []

for i in output:  # Output consists of angle value, thickness value, and reflection coefficient
    xcoord.append(i[0])
    ycoord.append(i[1])
    zcoord.append(i[2])


# Creating the array of analytical values expected
xcoord2 = []
ycoord2 = [] 
zcoord2 = []

for val in analytic:
    xcoord2.append(val[0])
    ycoord2.append(val[1])
    zcoord2.append(val[2])

ycoords = ycoord[::10]
xcoords = xcoord[::10]
zcoords = zcoord[::10]

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel("incidence angle (rad)")
ax.set_ylabel("thickness (nm)")
ax.set_zlabel("reflection")
#ax.set_title("Reflection spectrum for MgF2 layer on top of glass substrate")

ax.plot(xcoord2, ycoord2, zcoord2, color='k', zorder=1, label="analytical formula")
ax.scatter(xcoord, ycoord, zcoord, c=zcoord, cmap=cm.viridis, zorder=2, alpha=0.05)

ax.legend()

plt.show()
#  Checking if r + t = 1 (it should be equal to 1 as we haven't considered
# absorption yet).

print('Maximum t+r value = ' + str(max(tot_amp)) + ', minimum value =' + str(min(tot_amp)) + \
      ' Both of these values are on the order of 10^-16, suggesting this this a floating ' \
      'point error.')

plt.figure()
angle = 0.2

for j in d:  # numerically getting all values

    n1 = tmm.complx_n(incominglam, *MgF2)
    n2 = tmm.complx_n(incominglam, *BK7)
    n1 = complex(np.real(n1), 0)
    n2 = complex(np.real(n2), 0)
    ns = [n1, n2]
    ds = [j]

        # Obtaining the reflection and transmission for a range of values

    r, t, a = tmm.TMM(incominglam, angle, incomingpol, ns, ds, absorption=True)
    
    plt.scatter(j, r, color = 'k')

# Calculating the theoretical optimal parameters for each angle
ko = 2 * np.pi/incominglam
kx = ko * np.sin(angle)
        
d2 = np.pi/ (2 * np.sqrt((np.real(n1)*ko)**2 - kx**2))  # the analytical formula for d given the phase
        #print(np.shape(d2))
ds = [d2]

r2, t2 = tmm.TMM(incominglam, angle, incomingpol, ns, ds)
plt.scatter(d2, r2, label='Analytical optimal thickness')

plt.scatter(0, 0.040453, color='k', label='Datapoints computed using TMM')

plt.xlabel('Thickness of MgF2 layer (nm)')
plt.ylabel('Total reflection')
plt.title('Cross section of 3D plot for fixed angle '+str(angle))    
plt.legend()
plt.grid()


# %%
# Reflection and transmission in layer of gold.
# Investigating for normal incidience (angle = 0)

ang = 0
visible_spec = np.arange(500, 900, 1)
d_range = np.arange(10, 100, 1)
incomingpol = "s"

# Creating output list with: wavelength, thickness, reflection, and transmission
output_list = []

for lam in visible_spec:
    for d_val in d_range:
        n1 = tmm.complx_n(lam, *Au)  # refractive index of gold
        n2 = tmm.complx_n(lam, *BK7)  # refractive index of BK7 glass

        ns = [n1, n2]
        ds = [d_val]

        r, t, a = tmm.TMM(lam, ang, incomingpol, ns, ds, absorption=True)

        output_list.append((lam, d_val, r, t, a))

xcoord2, ycoord2, r_vals, t_vals, a_vals = [], [], [], [] ,[]

for output in output_list:
    xcoord2.append(output[0])
    ycoord2.append(output[1])
    r_vals.append(output[2])
    t_vals.append(output[3])
    a_vals.append(output[4])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel("Wavelegth")
ax.set_ylabel("thickness (nm)")
ax.set_zlabel("reflection")
ax.set_title("Reflection variation for Gold layer on top of glass substrate for normal incidence")

ax.scatter(xcoord2, ycoord2, r_vals, c=r_vals, cmap=cm.viridis, zorder=2, alpha=0.05)

ax.legend()

# Plotting the transmssion spectrum as a function of wavelength and thickness
fig2 = plt.figure()

ax2 = plt.axes(projection='3d')
ax2.set_xlabel("Wavelegth")
ax2.set_ylabel("thickness nm")
ax2.set_zlabel("Transmission coefficient")
ax2.set_title("Transmission spectrum for Gold layer on top of glass substrate for normal incidence")

ax2.scatter(xcoord2, ycoord2, t_vals, c=t_vals, cmap=cm.viridis, zorder=2, alpha=0.05)

ax2.legend()



# Transmission for fixed wavelength of 500nm as a function of layer thickness.

fixed_wavelength = 700
ang = 0  # normal incidence
data2d = []
for d_val in d_range:
    n1 = tmm.complx_n(fixed_wavelength, *Au)  # refractive index of gold
    n2 = tmm.complx_n(fixed_wavelength, *BK7)  # refractive index of BK7 glass

    ns = [n1, n2]
    ds = [d_val]
    r, t, a = tmm.TMM(fixed_wavelength, ang, incomingpol, ns, ds, absorption=True)

    data2d.append((d_val, r, t, a))

data2d = np.array(data2d)

plt.figure()
plt.scatter(data2d[:, 0], data2d[:, 2], label='Transmission spectrum for fixed wavelength '+str(fixed_wavelength)+' nm')
plt.title('Gold layer with glass substrate')
plt.xlabel('Thickness of layer (nm)')
plt.ylabel('Transmission coefficient')
plt.legend()
plt.grid()

# Plotting the observed absorption
plt.figure()
plt.scatter(data2d[:, 0], data2d[:, 3], label='Absorption spectrum for fixed wavelength '+str(fixed_wavelength)+' nm', color='k')
plt.title('Gold layer with glass substrate')
plt.xlabel('Thickness of layer (nm)')
plt.ylabel('Total power absorbed')
plt.legend()
plt.grid()


#%% 
#Testing thickness of DBR layers against reflectivity
Ta2O5 = np.loadtxt("Ta2O5.csv", skiprows=1, unpack=True, delimiter=",")
MgF2 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
Ta2O5[0] = Ta2O5[0] * 1000

visible_spec = np.arange(500, 900, 1)
polarisation = 'p'
incangle = 0
lam_opt = 633

N_range = np.arange(1,10,4)
test = np.arange(0, 2, 0.01)

plt.figure()
for N in N_range:
    r_values = []

    for i in test:
        dt = i*np.real(lam_opt / (2 * (tmm.complx_n(lam_opt, *Ta2O5))))
        dm = i*np.real(lam_opt / (2 * (tmm.complx_n(lam_opt, *MgF2))))
        n_stack2, d_stack2 = tmm.stacklayers(N, lam_opt, dm, dt, MgF2, Ta2O5)

        r, t = tmm.TMM(lam_opt, incangle, polarisation, n_stack2, d_stack2)
        r_values.append(r)

    plt.plot(test, r_values, label='N = ' + str(N))

plt.xlabel('Phase due to width of layer (*pi)')
plt.ylabel('Reflection coefficient ')
plt.title('Reflectivity spectrum as a function of phase for different N values')
plt.legend(loc='best')
plt.grid()

print("As expected, max reflectance for more than 1 period of DBR occurs at the 2*d*kz = pi ")
