import scipy as sp
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
    index2 = 1 + (a1*wavelength**2)/(wavelength**2 - b1) + (a2*wavelength**2)/(wavelength**2 - b2) +\
             (a3*wavelength**2)/(wavelength**2 - b3)

    index = sp.sqrt(index2)
    return index


wavelength, n, k= sp.loadtxt("BK7.txt", skiprows = 1, unpack = True)
print(wavelength, n, k )
wavelength = sp.array(wavelength) # to avoid problems with curvefit make everything into arrays
n = sp.array(n)
k = sp.array(k)

googled = list([1.03961212,0.231792344,1.01046945,6.00069867*10**3,2.00179144*10**4,1.03560653*10**8])

#popt, pcov = curve_fit(sellmeier, wavelength, n, p0 = [1,0,1,0,0,100]) # initial guess taken from internet

trial = sp.arange(330, 2500)
interplovals = list(map(lambda x: linterpol(x, wavelength, n), trial))
sellmeiervals = list(map(lambda x: sellmeier(x, *googled), trial))

plt.plot(trial, interplovals)
plt.plot(trial, sellmeiervals)
plt.show()

residuals = []

for i in range(len(trial)):
    residuals.append(interplovals[i] - sellmeiervals[i])

print("The averaage residual is " + str(sp.mean(residuals)))

plt.scatter(trial, residuals)
plt.show()

plt.plot(wavelength, n)
plt.show()

plt.plot(wavelength, k)
plt.show()