"""
Created on Wed Mar  3 00:43:27 2021

@author: leonardobossi1
"""
import numpy as np
import tmmfile as tmm
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.animation import FuncAnimation

# Set the default color cycle
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["r", "k", "c"])

# Importing data for refractive indices of different materials
MgF2 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
BK7 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
Au = np.loadtxt('Au.txt', skiprows = 1, unpack = True)
Ta2O5 = np.loadtxt("Ta2O5.csv", skiprows=1, unpack=True, delimiter=",")
TiO2 = np.loadtxt("Devore-o.csv", skiprows=1, unpack=True, delimiter=",")
SiO2 = np.loadtxt("Malitson.csv", skiprows=1, unpack=True, delimiter=",")

# Converting from micrometers to nanometers
Ta2O5[0] = Ta2O5[0] * 1000
TiO2[0] = TiO2[0] * 1000
SiO2[0] = SiO2[0] * 1000

# Adding a k column to TiO2 and SiO2 so they fit the expected format
TiO2 = list(TiO2)
SiO2 = list(SiO2)
Timgs = list(np.zeros(len(TiO2[0])))
Simgs = list(np.zeros(len(SiO2[0])))
TiO2.append(Timgs)
SiO2.append(Simgs)

fixed_wavelength = 633
incangle = 0
polarisation = "s"

# Finding the number of periods required to reach 99.99% reflectivity, d1=d2=50nm,
# We used materials Ta2O5 and MgF2 with air substrate.
visible_spec = np.arange(500, 900, 1)
lam_opt = 633

# Defining the optimal thicknesses for maximum reflection 
d1opt = np.real(lam_opt / (4 * (tmm.complx_n(lam_opt, *Ta2O5))))
d2opt = np.real(lam_opt / (4 * (tmm.complx_n(lam_opt, *MgF2))))

N, r_current, plot = tmm.find_N(0.9999, fixed_wavelength, d1opt, d2opt, incangle, polarisation, Ta2O5, MgF2)
nplot = []
rplot = []
aplot = []

for i in plot:
    nplot.append(i[0])
    rplot.append(i[1])
    aplot.append(i[2])
                   
plt.figure()
plt.xlabel("Number of stacks")
plt.ylabel("Reflectance and absorption ratio")
plt.title('Number of stacks (each stack has 2 layers) required to reach 99.99% reflectivity.')
plt.plot(nplot, rplot, label = "reflection", color='k')
plt.grid()
plt.legend()

# Plotting absorption
plt.figure()
plt.plot(nplot, aplot, label = "absorption")
plt.xlabel("Number of stacks")
plt.ylabel("Absorption ratio")
plt.title('Absorption in the stack preventing 99.99% reflectivity')
plt.grid()
plt.legend()

# Finding the reflectivity of gold at wavelength = 633 nm (task 12)
n_gold = tmm.complx_n(633, *Au)
n_substrate = tmm.complx_n(633, *BK7)
d_opt = [633 / (np.real(n_gold) * 4)]

r_gold = tmm.TMM(633, 0, 's', [n_gold, n_substrate], d_opt)

print('reflectivity of gold at 633nm wavelength = ' + str(r_gold[0]))

#%%
# Investigating different incoming angles

# Creating a stack
N_stack = 2  # Period of layers
d1 = 50
d2 = 50  # thicknesses of the two material layers

n_stack, d_stack = tmm.stacklayers(N_stack, 500, d1opt, d2opt, Ta2O5, MgF2)

ang_range = np.arange(0, np.pi / 2, 0.01)
r_output = []
t_output = []
angles = []

for ang in ang_range:
    r, t = tmm.TMM(fixed_wavelength, ang, polarisation, n_stack, d_stack)

    r_output.append(r)
    t_output.append(t)
    angles.append(ang)

# Plotting transmission and reflection
plt.figure()
plt.scatter(angles, r_output, label='R coefficient')
plt.xlabel("Angle of incidence (rad)")
plt.ylabel("Reflectance/Transmission coefficient")
plt.scatter(angles, t_output, label='T coefficient')
plt.title('R,T as a function of incident angle ' \
          ' for 4 layers of Ta2O5 and MgF2 with d=50nm')
plt.legend()
plt.grid()

# Setting N=14 and thicknesses as the optimal thickness for a wavelength of 633nm.

plt.figure()
angles = [0, np.pi / 6, np.pi / 3]

for ang in angles:

    r_values = []

    for lam in visible_spec:
        n_stack2, d_stack2 = tmm.stacklayers(14, lam, d1opt, d2opt, Ta2O5, MgF2)

        r, t = tmm.TMM(lam, ang, polarisation, n_stack2, d_stack2)
        r_values.append(r)

    if ang != 0:
        plt.plot(visible_spec, r_values, label='Incident Angle = \u03C0 /' + str(int(1 / ang * (np.pi))))
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
angle_tolim = [2 * np.pi / 6, 11 * np.pi / 24, 99 * np.pi / 200]

plt.figure()

for ang in angle_tolim:
    r_values = []
    for lam in visible_spec:
        n_stack2, d_stack2 = tmm.stacklayers(14, lam, d1opt, d2opt, Ta2O5, MgF2)

        r, t = tmm.TMM(lam, ang, polarisation, n_stack2, d_stack2)
        r_values.append(r)

    plt.plot(visible_spec, r_values, label='Incident Angle = \u03C0 /' + str(round(1 / ang * (np.pi),5)))

plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient ')
plt.title('Reflectivity as angle tends towards \u03C0 / 2 ')
plt.legend(loc='best')
plt.grid()

#%%
# Altering the N values to see how increasing layers affects the band gap.
# Note that the R spectrum as a function of wavelength,
# we use the optimal thicknesses for an incoming wavelength of 633nm
# so this is where the peak should be observed.

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color='kbgrcmy')
plt.figure()

N_range = np.arange(1, 10, 2)
incangle = 0


for N in N_range:
    r_values = []

    for lam in visible_spec:
        n_stack2, d_stack2 = tmm.stacklayers(N, lam, d1opt, d2opt, Ta2O5, MgF2)
        r, t = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2)
        r_values.append(r)

    plt.plot(visible_spec, r_values, label='N = ' + str(N))

plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflection coefficient ')
plt.title('Reflectivity spectrum as a function of wavelength for different N values')
plt.legend(loc='best')
plt.grid()

# Plotting spectrum from N=14.
plt.figure()

N = 14
r_values = []
a_values=[]
t_values=[]

layer1 = Ta2O5
layer2 = MgF2

# Ensuring the thicknesses are optimal for reflectivity
d1 = np.real(lam_opt / (4 * (tmm.complx_n(lam_opt, *layer1))))
d2 = np.real(lam_opt / (4 * (tmm.complx_n(lam_opt, *layer2))))

for lam in visible_spec:

    
    n_stack2, d_stack2 = tmm.stacklayers(N, lam, d1, d2, layer1, layer2)
    
    r, t, a = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2, absorption=True)
    r_values.append(r)
    a_values.append(a)
    t_values.append(t)

plt.plot(visible_spec, t_values, label='R')
#plt.plot(visible_spec, a_values, label='Absorption values')

# Creating a high/low edge pass filter by eliminating noise

N = 14
r_noise = []
t_noise = []

# Ensuring the thicknesses are optimal for reflectivity
d1 = np.real(lam_opt / (4 * (tmm.complx_n(lam_opt, *layer1))))
d2 = np.real(lam_opt / (4 * (tmm.complx_n(lam_opt, *layer2))))

for lam in visible_spec:
    

    n_stack2, d_stack2 = tmm.stacklayers(N, lam, d1, d2, layer1, layer2)
    d_stack2[0] = d_stack2[0]/2
    n_stack2[-1] = tmm.complx_n(lam, *layer1)
    n_stack2.append(1)
    d_stack2.append(d1/2)
    
    r, t, a = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2, absorption=True)
    r_noise.append(r)
    t_noise.append(t)


plt.plot(visible_spec, t_noise, label='R with noise reduction')
#plt.plot(visible_spec, a_values, label='Absorption values')

#plt.title('Noise reduction for optical filter applications')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission')
plt.legend(loc='lower right')
plt.grid()
plt.figure()

# Extending the stop-band
r_extend = [] 
t_extend = []

lam_opt2 = 736
d3 = np.real(lam_opt2 / (4 * (tmm.complx_n(lam_opt2, *layer1))))
d4 = np.real(lam_opt2 / (4 * (tmm.complx_n(lam_opt2, *layer2))))

for lam in visible_spec:
    
    n_1, d_1 = tmm.stacklayers(N, lam, d1, d2, layer1, layer2)
    n_2, d_2 = tmm.stacklayers(N, lam, d3, d4, layer1, layer2)
    
    d_1[0] = d_1[0]/2
    d_2[0] = d_2[0]/2
    n_2[-1] = tmm.complx_n(lam, *layer1)
    n_1[-1] = tmm.complx_n(lam, *layer1)
    n_2.append(1)
    d_2.append(d3/2)
    d_1.append(d1/2)
    
    n_stack2 = np.append(n_1, n_2)
    d_stack2 = np.append(d_1, d_2)
    
    
    r, t, a = tmm.TMM(lam, incangle, polarisation, n_stack2, d_stack2, absorption=True)
    r_extend.append(r)
    t_extend.append(t)

plt.plot(visible_spec, t_extend, label='Extended stopband', color='r')
plt.plot(visible_spec, t_noise, label='Noise reduction')

plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission')
#plt.title('Tranmission spectrum after extending the stop band and reducing noise')
plt.legend(loc='upper left')
plt.grid()



#%%
# Investigating changing materials and the effect it has on spectral width

# Using a stack with different materials (SiO2 with MgF2, TiO2 with Ta2O5)
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["r", "k", "c"])

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

N = 60
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
# plt.savefig("High resoltion.png",dpi=300)

print('spectral width for SiO2 and MgF2 with N=60 is ' + str(tmm.spectral_width(r_values3)) + ' nm')
print('spectral width for TiO2 with SiO2 with N=14 is ' + str(tmm.spectral_width(r_values2)) + ' nm')
print('spectral width for TiO2 and Ta2O5 with N=14 is ' + str(tmm.spectral_width(r_values)) + ' nm')

# %%
# for wavelength of 633nm, phase response of reflected wave with both gold and dbr
n_gold = tmm.complx_n(633, *Au)
n_substrate = tmm.complx_n(633, *BK7)
dG = [633 / (np.real(n_gold) * 4)]

nT = tmm.complx_n(fixed_wavelength, *Ta2O5)
nM = tmm.complx_n(fixed_wavelength, *MgF2)
dT = 633 / (np.real(nT) * 4 * np.cos(incangle))
dM = 633 / (np.real(nM) * 4 * np.cos(incangle))

n_stack, d_stack = tmm.stacklayers(N, 633, dM, dT, MgF2, Ta2O5, substrate_n=n_substrate)

var_wavelength = np.arange(553, 741)
rsphases = []
rgphases = []

for i in var_wavelength:
    n_stack, d_stack = tmm.stacklayers(14, i, dM, dT, MgF2, Ta2O5, substrate_n=n_substrate)
    rstack, tstack = tmm.TMM(i, 0, "s", n_stack, d_stack, squared=False)
    rgold, tgold = tmm.TMM(i, 0, "s", [n_gold, n_substrate], dG, squared=False)
    rsphase = np.angle(rstack)
    rgphase = np.angle(rgold)
    rsphases.append(rsphase)
    rgphases.append(rgphase)

plt.figure()
plt.plot(var_wavelength, rsphases, label="Phase response of DBR")
plt.plot(var_wavelength, rgphases, label="Phase response of gold")
plt.grid()

plt.legend()
plt.xlabel("Wavelength in nm")
plt.ylabel("Phase of reflected wave in rad")
plt.yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
           [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$'])

# %%
# Introducing a defect to see what happens
# expansion of a central cavity
# for DBR of Ta2O5, MgF2
# central optical cavity of air to begin
nT = tmm.complx_n(fixed_wavelength, *Ta2O5)
nM = tmm.complx_n(fixed_wavelength, *MgF2)
dT = 633 / (np.real(nT) * 4 * np.cos(incangle))
dM = 633 / (np.real(nM) * 4 * np.cos(incangle))

ncav = complex(1 + 0j)
dcavs = np.arange(0, 500, 10)

n_stack, d_stack = tmm.stacklayers(15, 633, dM, dT, MgF2, Ta2O5, substrate_n=n_substrate)

var_wavelength = np.arange(500, 900, 0.05)
r_output = []
animdata = []
for dcav in dcavs:
    r_output = []
    for i in var_wavelength:
        n_stack, d_stack = tmm.stacklayers(14, i, dM, dT, MgF2, Ta2O5, substrate_n=n_substrate)
        n_stack1 = n_stack.copy()
        n_stack1.pop()
        n_stack1.append(ncav)
        n_stack1.extend(n_stack)
        d_stack1 = d_stack.copy()
        d_stack1.append(dcav)
        d_stack1.extend(d_stack)
        rstack, tstack = tmm.TMM(i, 0, "s", n_stack1, d_stack1)
        r_output.append(rstack)
    print(dcav)
    #    im = plt.figure()
    #   plt.plot(var_wavelength, r_output)
    #  plt.title("Cavity of air at " + str(dcav) + "nm thickness")
    # plt.xlabel('wavelength (nm)')
    # plt.ylabel("Reflectance")
    animdata.append(r_output)

fig, ax = plt.subplots(figsize=(5, 3))
ax.set(xlim=(500, 900), ylim=(0, 1.1))
line = ax.plot(var_wavelength, animdata[0], color='k', lw=2)[0]


def animate(i):
    line.set_ydata(animdata[i])


anim = FuncAnimation(
    fig, animate, interval=100, frames=len(dcavs) - 1)

plt.draw()
plt.show()
anim.save('filename.gif', writer='imagemagick')

#%%
# Trying to make a band gap
# expansion of a central cavity
# for DBRs of Ta2O5, MgF2

nL = tmm.complx_n(fixed_wavelength, *MgF2)
dL = 633 / (np.real(nL) * 4 * np.cos(incangle))

nH = tmm.complx_n(fixed_wavelength, *Ta2O5)
dH = 633 / (np.real(nH) * 4 * np.cos(incangle))

ncav = 1.0003
dcavs = np.arange(0, 500, 0.5)

var_wavelength = np.arange(500, 900, 1)
r_output = []
t_output = []
animdata = []
animdatat = []
rout2 = []
drec = []
for i in var_wavelength:
    #n_cav = tmm.complx_n(i, *GaAs)
    n2stack, d2stack = tmm.stacklayers(4, i, dL, dH, MgF2, Ta2O5, substrate_n=n_substrate)
    n1stack, d1stack = tmm.stacklayers(4, i, dH, dL, Ta2O5, MgF2, substrate_n=n_substrate)
    n1stack.pop()
    n2stack.pop()
    n_stack = n1stack.copy()
    n_stack.extend(n2stack)
    d_stack = d1stack.copy()
    d_stack.extend(d2stack)

    n_stack.append(tmm.complx_n(i, *MgF2))
    d_stack.append(dL)

    n_stack.extend(n1stack)
    n_stack.extend(n2stack)
    d_stack.extend(d1stack)
    d_stack.extend(d2stack)

    n_stack.append(tmm.complx_n(i, *BK7))

    rstack, tstack = tmm.TMM(i, 0, "s", n_stack, d_stack)
    t_output.append(tstack)
    r_output.append(rstack)


plt.figure()
plt.plot(var_wavelength, t_output)
plt.show()
