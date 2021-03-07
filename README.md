Read me for Transfer-Matrix-Method project Code.

Python version used: 3.7.6
Libraries used in code: Matplotlib, Numpy, mpl_toolkits.mplot3d, cmath.

This file includes the modules used in the investigation as well as the datasets for the complex refractive indices of different materials.
The complex refractive index data is obtained from https://refractiveindex.info and are either in .txt or .csv format.

The main modules include:

tmmfile.py
Inlcudes all the functions used in the analysis of the transmission/refractive spectrum. Importantly, it includes the code for the transfer matrix method.

test.py 
Testing the transfer matrix works as predicted using test cases we understand.

SPP.py
Investigating the behaviour of a gold layer on a glass substrate, including the surface plasmon polarisation and absorption spectrum of gold.

DBR.py
Investigating distributed bragg reflectors (DBR) and using them to create optical filters. 

For all of the modules, tmmfile.py is imported as the functions in the file are needed for analysis.




