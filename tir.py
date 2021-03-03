"""
Created on Wed Mar  3 00:43:27 2021

@author: leonardobossi1
"""



import numpy as np
import tmmfile as tmm
import matplotlib.pyplot as plt
import matplotlib as mpl

# Importing different materials
MgF2 = np.loadtxt("MgF2.txt", skiprows=1, unpack=True)
BK7 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)
Au = np.loadtxt('Au.txt', skiprows = 1, unpack = True)
Ta2O5 = np.loadtxt("Ta2O5.csv", skiprows=1, unpack=True, delimiter=",")

material = MgF2

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["b", "k", "r"])

# Finding total internal reflection
inc_lam = 633

n = tmm.complx_n(inc_lam, *material)

nglass = tmm.complx_n(inc_lam, *BK7)

critical_ang = np.arcsin(1.0003/np.real(n))

angles = np.arange(0, np.pi / 2, 0.01)

r_p = []
r_s = []

for ang in angles:
    rp, tp = tmm.TMM(inc_lam, ang, 'p', [1.0003], [], glass_inc=True, material = material)
    r_p.append(rp)

    rs, ts = tmm.TMM(inc_lam, ang, 's', [1.0003], [], glass_inc=True, material= material)
    r_s.append(rs)

plt.figure()
plt.plot(angles, r_p, label='p polarisation')
plt.plot(angles, r_s, label='s polarisation')
plt.xlabel('Angle (rad)')
plt.ylabel('Reflectance')
plt.title("Total Internal Reflection from MgF2 to air")

plt.vlines(critical_ang, 0, 1, color='r', label='critical angle')

plt.grid()
plt.legend()

#%%

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["b", "k", "r"])
# comparison between p (TM) polarisation and s
inc_lam = 633

d_gold = [50]
critical_ang = np.arcsin(1 / np.real(tmm.complx_n(inc_lam, *BK7)))

n_list = [tmm.complx_n(inc_lam, *Au), 1.003]

angles = np.arange(0, np.pi / 2, 0.001)

r_p = []
r_s = []

for ang in angles:
    rp, tp = tmm.TMM(inc_lam, ang, 'p', n_list, d_gold, glass_inc=True, material = BK7)
    r_p.append(rp)

    rs, ts = tmm.TMM(inc_lam, ang, 's', n_list, d_gold, glass_inc=True, material = BK7)
    r_s.append(rs)

plt.figure()
plt.plot(angles, r_p, label='p polarisation')
plt.plot(angles, r_s, label='s polarisation')
plt.xlabel('Angle (rad)')
plt.ylabel('Reflectance')
plt.title("Finding the Surface Plasmon Polariton")
plt.vlines(critical_ang, 0, 1, color='r', label='critical angle')

plt.grid()
plt.legend()

tmm.TMM(633, 0, 'p', n_list, [50], glass_inc=True,material =  BK7)

