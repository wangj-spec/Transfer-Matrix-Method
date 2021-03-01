

import numpy as np
import tmmfile as tmm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
plt.style.use('ggplot')

# Importing different materials
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

# Set the default color cycle
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["b", "k", "r"])

# comparison between p (TM) polarisation and s
inc_lam = 633
d_gold = [49]
critical_ang = np.arcsin(1.0003 / np.real(tmm.complx_n(inc_lam, *BK7)))

n_list = [tmm.complx_n(inc_lam, *Au), 1.0003]

angles = np.arange(0, np.pi / 2, 0.01)

r_p = []
r_s = []

for ang in angles:
    rp, tp = tmm.TMM(inc_lam, ang, 'p', n_list, d_gold, glass_inc=True, material=BK7)
    r_p.append(rp)

    rs, ts = tmm.TMM(inc_lam, ang, 's', n_list, d_gold, glass_inc=True, material=BK7)
    r_s.append(rs)

plt.figure()
plt.plot(angles, r_p, label='p polarisation')
plt.plot(angles, r_s, label='s polarisation')

plt.vlines(critical_ang, 0, 1, color='r', label='critical angle')

plt.grid()
plt.legend()
