#!/usr/bin/python

import numpy as np
from numpy import log,exp
#from scipy import *
import scipy.optimize as opt
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import axis as ax
import matplotlib               
from pylab import *
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

x_t ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

eiz_x = np.loadtxt('data/eiz_x_output.txt') #parallel
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') #perpendicular

# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV #(2.41*1e14) # in rad/s
ns = arange(0.0,500.0)
z = ns * coeff
kb = 8.6173e-5 #in eV/K
T = 300.0 # r.t. in K
c = 3.0e17
r_1 = 0.5
r_2 = 0.5
t = np.linspace(0.0,1000.0,100)
G_const = -1.0*r_1*r_2
#ls = 8e-10#[1e-6 , 1e-5 , 1e-4]
ls = linspace(0.1 , 4.0 , 50)
thetas = linspace(np.pi/4.0,(7./8)*np.pi,6)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
eiz_m = 1.0

height= np.zeros(shape=(len(ls),len(thetas)))
aiz = np.zeros(len(eiz_z))
def a(perp, par):
    return 2.0*(perp-1.0)/((perp+1.0)*(par-1.0))
aiz = a(eiz_z,eiz_x)

#
lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']

int_g_arg = np.empty(len(ns))
int_g = np.empty(len(ns))
A = np.empty(len(ns))
A0 = np.empty(len(ns))
A2 = np.empty(len(ns))
sum_A = np.empty(len(ls))
G = np.empty(len(ls))

pl.figure()
for p in range(len(thetas)):
	for k in range(len(ls)):
		for j in range(len(ns)):
			int_g[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(((1.0 - aiz[j])*(1.0 - aiz[j])*(t*t + 2.0)*(t*t + 2.0))),t)

			A[j] = 1.0/32.0*((eiz_x[j]-1.0)/1.0)*\
			((eiz_x[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g[j]	
			
			sum_A[k] = sum(A[j])
		print sum_A	
	G = -((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A

	savetxt("data/G2_output.txt", G)
	pl.plot(ls,G, color = colors[p])#, linestyle = lts[k])	
#figure()
#contourf(X,Y,height,100); clim(-1e-47,0)#,cmap = hot())
#levels = [0.1,0.5,1.0,1.2,1.3,1.4,1.5,2.0]#,0.80,1.00]
######contour(X,Y,height)
#CS = contour(X,log(Y),height)#,levels, colors = 'k')
#man_loc = [(1,1),(2,2),(3,3),(4,4)]
#clabel(CS, inline =1,fmt = '%2.1f', fontsize = 18,color = 'k', manual = man_loc)
#contourf(X,Y,height, 1000); clim(0,4)#;cbar = colorbar();
#CS=contour(X,Y,height, levels)#, colors = 'k')
#clabel(CS,fmt = '%2.1f',colors = 'k',fontsize = 18)
#xlabel(r'$\mathrm{\Phi_{PEG\,400}\,\,added\, to\, fixed\,\,\Phi_{PEG\,3500}=\,0.20}$', size = 18)#\,=\,\frac{2}{\pi}C$', size = 20)
#ylabel(r'$\Delta f$', size = '20')#\,\,\,Energy cost to enter pore$', size = 20)
#cbar.ax.set_ylabel(r'$P\,=\,\frac{\Phi_{400}^{in}}{\Phi_{400}^{out}}$', size = 20)
#cbar.add_lines(CS)
#axes([0,0.78,0,3])
pl.show()
#####savefig('plots/130930_test')
#cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
#cbar.add_lines(CS)
#####show()
####close()



