#!/usr/bin/python
import matplotlib               
#import pyreport
import numpy as np                  
from pylab import *
#from pylab import show
from matplotlib import pyplot as pl

#x_x,y_x = np.loadtxt('data/CNT6_5_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#x_z,y_z = np.loadtxt('data/CNT6_5_ze2_solid_30.txt',unpack=True, usecols = [0,1])

#x_x,y_x = np.loadtxt('data/CNT9_0_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#x_z,y_z = np.loadtxt('data/CNT9_0_ze2_solid_30.txt',unpack=True, usecols = [0,1])

#x_x,y_x = np.loadtxt('data/CNT9_1_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#x_z,y_z = np.loadtxt('data/CNT9_1_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#
#x_x,y_x = np.loadtxt('data/CNT9_3_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#x_z,y_z = np.loadtxt('data/CNT9_3_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#
x_x,y_x = np.loadtxt('data/CNT29_0_xe2_solid_30.txt',unpack=True, usecols = [0,1])
x_z,y_z = np.loadtxt('data/CNT29_0_ze2_solid_30.txt',unpack=True, usecols = [0,1])

x_w,y_w = np.loadtxt('data/water-L.txt',unpack=True, usecols = [0,1])

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#------------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV #(2.41*1e14) # in rad/s
T = 300.0
n = arange(0,500)
z = n * coeff

eiz_x = empty(len(z))
eiz_z = empty(len(z))
eiz_w = empty(len(z))

eiz_x_arg=empty(len(x_x))
eiz_z_arg=empty(len(x_z))
eiz_w_arg=empty(len(x_w))

for j in range(len(z)):
    for i in range(len(x_x)):
        eiz_x_arg[i]=x_x[i]*y_x[i] / (x_x[i]**2 + z[j]**2)
    eiz_x[j] = 1 + (2./pi) * trapz(eiz_x_arg,x_x)

    for m in range(len(x_z)):
        eiz_z_arg[m]=x_z[m]*y_z[m] / (x_z[m]**2 + z[j]**2)
    eiz_z[j] = 1 + (2./pi) * trapz(eiz_z_arg,x_z)    

    for p in range(len(x_w)):
        eiz_w_arg[p]=x_w[p]*y_w[p] / (x_w[p]**2 + z[j]**2)
    eiz_w[j] = 1 + (2./pi) * trapz(eiz_w_arg,x_w)    

#savetxt("data/eiz_x_65.txt", eiz_x)
#savetxt("data/eiz_z_65.txt", eiz_z)

#savetxt("data/eiz_x_90.txt", eiz_x)
#savetxt("data/eiz_z_90.txt", eiz_z)

#savetxt("data/eiz_x_91.txt", eiz_x)
#savetxt("data/eiz_z_91.txt", eiz_z)
#
#savetxt("data/eiz_x_93.txt", eiz_x)
#savetxt("data/eiz_z_93.txt", eiz_z)
#
savetxt("data/eiz_x_290.txt", eiz_x)
savetxt("data/eiz_z_290.txt", eiz_z)
#
#savetxt("data/eiz_w.txt", eiz_w)

#pl.figure()
#pl.plot(x_t,y_t,    color = 'k', label = 'total')
#pl.plot(x_x+10,y_x, color = 'b', label = r'$\hat{x}$')
#pl.plot(x_y+20,y_y, color = 'g', label = r'$\hat{y}$')
#pl.plot(x_z+30,y_z, color = 'r', label = r'$\hat{z}$')
#pl.xlabel(r'$\hbar\omega\,\,\,[eV]\,\,\,!shifted\,10\,eV\,for\,visualization$', size = 24)
#pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
#pl.legend()
#pl.title(r'CG-10 DNA eps2... shifted for visualization')
##pl.savefig('plots/131010_Hopkins_CG10_eps2_x_z.png', dpi = 300 )
#pl.savefig('plots/131010_Hopkins_CG10_eps2_all_directions.pdf')
##pl.show()
#
pl.figure()
#pl.plot(x_t,y_t, color = 'k', label = r'$\varepsilon^{\prime\prime}_{total}(\omega)$')
pl.plot(x_x,y_x, color = 'b', label = r'$\varepsilon^{\prime\prime}_\hat{x}(\omega)$')
#pl.plot(x_y,y_y, color = 'g', label = r'$\varepsilon^{\prime\prime}_\hat{y}(\omega)$')
pl.plot(x_z,y_z, color = 'r', label = r'$\varepsilon^{\prime\prime}_\hat{z}(\omega)$')
pl.plot(x_w,y_w, color = 'c', label = r'$\varepsilon^{\prime\prime}_{H_{2}O}(\omega)$')
pl.xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
pl.legend()
#pl.title(r'[6,5] and water eps2')
#pl.savefig('plots/140306_65w65_eps2.pdf')
pl.show()
#pl.close()
#imshow(1)
#pl.savefig('plots/131010_Hopkins_CG10_eps2_x_z.png', dpi = 300 )
#
##pl.figure()
###pl.plot(x_t,y_t, label = 'total')
###pl.plot(x_x,y_x, label = r'$\hat{x}$')
##pl.plot(x_y,y_y, label = r'$\hat{y}$')
##pl.plot(x_z,y_z, label = r'$\hat{z}$')
##pl.xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
##pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
##pl.legend()
##pl.show()
###pl.close()
###imshow(1)
##pl.savefig('plots/DNA_spectra_x_y.png', dpi = 300 )
#
pl.figure()
pl.plot(n,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
#pl.plot(n,eiz_y, color = 'g', label = r'$\varepsilon_{\hat{y}}(i\zeta_{N})$')
pl.plot(n,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{n})$')
pl.plot(n,eiz_w, color = 'c', label = r'$\varepsilon_{\hat{w}}(i\zeta_{n})$')
pl.xlabel(r'$N$', size = 24)
pl.ylabel(r'$\varepsilon(i\zeta)$', size = 24)
pl.legend()
#pl.title(r'[6,5] and water eiz')
#pl.savefig('plots/DNA_eiz_x_z.png', dpi = 300 )
#pl.savefig('plots/140306_65w65_eiz.pdf')
pl.show()
#pl.close()


