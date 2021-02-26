#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:28:54 2021

@author: leonardobossi1
"""

import numpy as np
import cmath as cm

BK7 = np.loadtxt("BK7.txt", skiprows=1, unpack=True)

def linterpol(wavelength, wavedata, indexdata):
    '''
    Parameters:
        wavelength:: float
        wavedata:: numpy.array
            Discrete wavelength values used to interpolate
        indexdata:: numpy.array
            Corresponding refractive index values
    Returns:
        Interpol::float
            Interpolated refractive index value for given wavelength.
    '''
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
    '''
    Parameters:
        lam:: float
            wavelength
        lam_data:: np.array
            array of wavelength values obtained from online resources
        real_data:: np.array
            corresponding real values of refractive index
        img_data:: np.array
            corresponding imaginary values of refractive index
    Returns:
        n:: complex
            complex refractive index at provided wavelength.
    '''
    
    if lam > lam_data[-1] or lam < lam_data[0]:
        raise Exception('Inputted value is out of range provided in the data')
        
    n_real = linterpol(lam, lam_data, real_data)
    n_img = linterpol(lam, lam_data, img_data)
    n = complex(n_real, n_img)
    
    return n


def chifunction (k0, kx, polarisation, ncomplex1, ncomplex2):
    '''
    Parameters:
        k0:: float
            Wave number in free space.
        kx:: float
            k in the x direction.
        polarisation:: string
            polarisation of wave (can be p or s polarised).
        ncomplex1:: complex
        ncomplex2:: complex
            Complex refractive indices of the two materials.
    Returns:
        chim:: float
            Chi_minus (to be used in the interfacial matrix in TMM.)
        chip:: float
            Chi_plus (to be used in the interfacial matrix in TMM.)

    '''

    kz1 = cm.sqrt((ncomplex1 * k0)**2 - kx**2)
    kz2 = cm.sqrt((ncomplex2 * k0)**2 - kx**2)
    
    if np.imag(kz2) < 0:
        kz2 = complex(np.real(kz2), -np.imag(kz2))
    if np.imag(kz1) < 0:
        kz1 = complex(np.real(kz1), -np.imag(kz1))
        

    if polarisation == "p":

        alpha = (ncomplex2/ncomplex1) * cm.sqrt(kz1/kz2)
        chip = (alpha + (1/alpha))/2
        chim = (alpha - (1/alpha))/2

    elif polarisation == "s":

        alpha = cm.sqrt(kz2/kz1)
        chip = (alpha + (1/alpha))/2
        chim = (alpha - (1/alpha))/2

    else:
        raise Exception("That is not a supported polarisation of light")
    
    return chip, chim  
    


def TMM(wavelength, angle, polarisation, ns, ds, squared = True, glass_inc = False):
    '''
    Parameters:
        wavelength:: float
            wavelength of incoming wave
        angle:: float
            incident angle.
        polarisation:: string
        ns:: list
            list of refractive indices for all the layers in the multi-layered stack.
        ds:: list
            list of all separations between the layers.
    Returns:
        r:: float
            Reflection coefficient.
        t:: float
            Transmission coefficient. 
    '''
    if polarisation != 's' and polarisation != 'p':
        raise Exception("That is not a supported polarisation of light")
        
    k0 = 2 * np.pi / wavelength
    if glass_inc:
        kx = complx_n(wavelength, *BK7) * k0 * np.sin(angle)
    else:
        kx = k0 * np.sin(angle)
        
    M = np.array([[1,0],[0,1]]) # initialise general matrix
    
    for i in range(len(ds)): # ds should be one item shorter than ns, the final ns should be for the substrate
        
        if i == 0: # setting incident medium to be air
            if glass_inc:
                n1 = complx_n(wavelength, *BK7)
            else:
                n1 = 1.0003 # air
            n2 = ns[i] # the refractive index of the layer

        else:
            n1 = ns[i-1]
            n2 = ns[i]

        kz = cm.sqrt((n2 ** 2 * k0 ** 2) - kx ** 2)
        
        if np.imag(kz) < 0:
            kz = complex(np.real(kz), -np.imag(kz))

        P_i = np.array([[np.exp(1j * kz * ds[i]), 0],
                        [0, np.exp(-1j * kz * ds[i])]]) # complex conjugate for backward propagation 
        
        chip, chim = chifunction(k0, kx, polarisation, n1, n2)
        T_i = np.array([[chip , chim],[chim, chip]])

        #interstep = np.matmul(P_i, T_i)
        M = np.matmul(P_i, np.matmul(T_i, M))

      
    n1 = ns[-2]
    n2 = ns[-1]

    chip, chim = chifunction(k0, kx, polarisation, n1, n2)
    T_i = np.array([[chip , chim],[chim, chip]]) # interfacial for the substrate
    M = np.matmul(T_i, M)

    r = -(M[1][0] / M[1][1])
    t = M[0][0] + M[0][1] * r

    r2 = abs(r) ** 2
    t2 = abs(t) ** 2  # want the transmittance and reflectance

    if squared == True:  # returns transmittance and reflectance for power
        return r2, t2
    
    elif squared == False:
        return r, t  # returns fresnel coefficients of r, t for electric field


def stacklayers(N, wavelength, d1 , d2, material1, material2, substrate_n = 1):
    ns = []
    ds = []
    n1 = complx_n(wavelength,*material1)
    n2 = complx_n(wavelength,*material2)


    for i in range(N):
        ns.append(n1)
        ds.append(d1)
        ns.append(n2)
        ds.append(d2)
        
    ns.append(substrate_n) 

    return ns, ds


def find_N(r_val, wavelength, d1, d2, angle, polarisation, material1, material2):
    N = 1
    plot = []
    r_current = 0 # initialies r_current
    
    while N < 15:
        
        ns, ds = stacklayers(N, wavelength, d1, d2, material1, material2)
        r_current = TMM(wavelength, angle, polarisation, ns, ds, squared=True)[0] 
        print(r_current) 
        plot.append([N, r_current])
        N += 1
        
    plot.append([N, r_current])

    return N, r_current, plot


