#!/usr/bin/python
import matplotlib               
#import pyreport
import numpy as np                  
from pylab import *
#from pylab import show
from matplotlib import pyplot as pl

x_x,y_x_unsc = np.loadtxt('data/CNT29_0_xe2_solid_30.txt',unpack=True, usecols = [0,1])
x_z,y_z_unsc = np.loadtxt('data/CNT29_0_ze2_solid_30.txt',unpack=True, usecols = [0,1])
x_w,y_w      = np.loadtxt('data/water-L.txt',unpack=True, usecols = [0,1])

y_x = y_x_unsc*1.#4.949
y_z = y_z_unsc*1.#4.949

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#------------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV #(2.41*1e14) # in rad/s
#coeff = 2.41e14 # in (1 rad)*(1/s)=inverse seconds
T = 297.0
#kb_J = 1.3806488e-23 # in J/K
#hbar = 6.625e-34 # in J/s
#coeff_J = 2.0*np.pi*kb_J*T/hbar#1.602e-19*0.159e15 # in eV #(2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff
#coeff_J = 1.602e-19*0.159e15 # in eV #(2.41*1e14) # in rad/s

#z = n * coeff
#z = n * coeff_J

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
#
savetxt("data/eiz_x_output_eV.txt", eiz_x)
savetxt("data/eiz_z_output_eV.txt", eiz_z)
savetxt("data/eiz_w_output_eV.txt", eiz_w)

a =  Aiz(eiz_x,eiz_z,eiz_w)

pl.figure()
pl.plot(x_x,y_x, color = 'b', label = r'$\varepsilon^{\prime\prime}_\hat{x}(\omega)$')
pl.plot(x_z,y_z, color = 'r', label = r'$\varepsilon^{\prime\prime}_\hat{z}(\omega)$')
#pl.plot(x_z,y_z, color = 'r', label = r'$first\,peak:\,\,%6.2f$'%max(y_z))
pl.plot(x_w,y_w, color = 'c', label = r'$\varepsilon^{\prime\prime}_{H_{2}O}(\omega)$')
pl.axis([0,35,0,25])
pl.xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
pl.legend()
pl.title(r'[29,0] and water eps2')
pl.savefig('plots/290w290_eps2.pdf')
pl.show()
#
fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(n,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
ax.plot(n,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{n})$')
ax.plot(n,eiz_z, color = 'r', label = r'$max\,%6.2f$'%max(eiz_z))
ax.plot(n,eiz_w, color = 'c', label = r'$\varepsilon_{\hat{w}}(i\zeta_{n})$')
pl.axis([0,500,0,10])
pl.xlabel(r'$N$', size = 24)
pl.ylabel(r'$\varepsilon(i\zeta)$', size = 24)
#pl.legend()
pl.title(r'[29,0] and water eiz')

ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
ax_inset.plot(n, a,'k')#, linewidth = 2)#,label=r'$a(i\xi_{N})$')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$N$', size = 14)
pl.ylabel(r'$a(i\xi)$', size = 14)
pl.savefig('plots/290w290_eiz.pdf')
pl.show()

