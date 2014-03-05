#!/usr/bin/python
import matplotlib               
#import pyreport
import numpy as np                  
from pylab import *
#from pylab import show
from matplotlib import pyplot as pl

x_t ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

eiz_x = np.loadtxt('data/eiz_x_output.txt',unpack=True, usecols = [0])
eiz_z = np.loadtxt('data/eiz_z_output.txt',unpack=True, usecols = [0])

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#------------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV #(2.41*1e14) # in rad/s

n = arange(0,250)
z = n * coeff
diff_eiz = eiz_z - eiz_x
diff_eps = y_z - y_x
#################################################################

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(n, eiz_x, linewidth = 2)#, label = r'Synthesized')#$ \epsilon (i \zeta ) _{peak}$')
ax.plot(n, eiz_z, linewidth = 2)#, label = r'Ab initio')#$ \epsilon (i \zeta )_{no \,peak} $')
pl.xlabel(r'$N$', size = 24)
pl.ylabel(r'$ \epsilon (i \zeta_{N} ) $', size = 24)
#pl.axis([0.0, 150,1.3,2.4])

ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
ax_inset.plot(n,diff_eiz, linewidth = 2,color= 'r',linestyle='-')
#pl.title('Difference $\epsilon$(i$\zeta$)',size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$N$',size = 14)#(r'$\zeta$')
pl.ylabel(r'$ \delta \epsilon (i  \zeta_{N} )$',size = 14)
#pl.axis([0.0, 150,0.0,0.1])
#pp.savefig()
pl.savefig('plots/131006_delta_eiz.png', dpi = 300)
#pl.savefig('130807_plots/fig2b_eiz_deiz.pdf')
#####################################################################
pl.figure()
pl.plot(x_x,y_x, linewidth = 2)#, label = r'$\hat{x}$')
pl.plot(x_z,y_z, linewidth = 2)#, label = r'$\hat{z}$')
pl.xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
pl.legend()
pl.savefig('plots/131006_eps2_x_z.png', dpi = 300)

pl.figure()
pl.plot(x_x,diff_eps, linewidth = 2)#, label = r'$\hat{x}$')
pl.xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
pl.legend()
pl.show()

