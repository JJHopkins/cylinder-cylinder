#!/usr/bin/python
import matplotlib               
import pyreport
import numpy as np                  
from pylab import *
from pylab import show
from matplotlib import pyplot as pl
#$
#$ $y = x^2$:
#$
x_t ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

z_t = np.ma.polyfit(x_t, y_t, 3)#deg, rcond=None, full=False, w=None, cov=False)
#z_t = np.polynomial.polynomial.polyfit(x_t, y_t,10)#deg, rcond=None, full=False, w=None, cov=False)
z_t_fit = np.poly1d(z_t)
x_fit = np.linspace(0,50,50)

z_z = np.ma.polyfit(x_z, y_z, 5)#deg, rcond=None, full=False, w=None, cov=False)
z_z_fit = np.poly1d(z_z)
x_fit = np.linspace(0,50,25)
y_interp = np.interp(x_fit, x_t,y_t)
#! trying to see if I can show this plot in pyreport
pl.figure()
pl.plot(x_t,y_t, label = 'total')
pl.plot(x_x,y_x, label = r'$\hat{x}$')
pl.plot(x_y,y_y, label = r'$\hat{y}$')
pl.plot(x_z,y_z, label = r'$\hat{z}$')
pl.plot(x_fit,z_t_fit(x_fit), label = r'$t\, fit$')
pl.plot(x_fit,y_interp, label = r'$t\, interp$')
pl.xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
pl.legend()
show()
#pl.close()
#imshow(1)
#pl.savefig('plots/DNA_spectr_tot_x_y_z.pdf')
