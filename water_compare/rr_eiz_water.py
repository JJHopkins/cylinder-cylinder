
#!/usr/bin/env python
import numpy as np
from scipy.integrate import romb
from scipy.integrate import simps #quad
from scipy.integrate import quad
import matplotlib.pyplot as pl

# Input dielectric response data
x,eiz_w_rr = np.loadtxt('LD_water_pars1.txt', delimiter='     ', usecols=(0,1), unpack=True)

coeff = 2.411e14 # in rad/s
ns = np.arange(0.,500.)
z = ns * coeff
# NOTES:
T = 300
# h_bar = 1. #1.0546e-34 #in Js
kb = 8.6173e-5  # in eV/K
# at RT, 1 kT = 4.11e-21 J
# 1 eV = 1.602e-19 J = 0.016 zJ
h_bar = 6.5821e-16 #eVs
z_eV = (2*np.pi*kb*T/h_bar)*ns 
#	= (0.159 eV) / (6.5821e-16 eVs) 
#	= n*2.411e14 rad/s
# z_n_J = (2*pi*kT/h_bar)n 
#	= (1.3807e-23 J/K) / (1.0546e-34 Js))*n
#	= n*2.411e14 rad/s
#coeff = 0.159 # in eV w/o 1/h_bar

#a = range(0,16001)

pl.figure()
pl.plot(x,eiz_w_rr)
#pl.plot(ns,z)
#pl.plot(ns,z_eV)

pl.show()
