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

x_t_eV ,y_t  = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x_eV ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y_eV ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z_eV ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

x_t =x_t_eV /6.582e-16
x_x =x_x_eV /6.582e-16
x_y =x_y_eV /6.582e-16
x_z =x_z_eV /6.582e-16

eiz_x = np.loadtxt('data/eiz_x_output.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') # parallel,axial

# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 2.41e14 # in rad/s
ns = arange(0.0,500.0)
T = 300.0 # r.t. in K
kb_J =1.3806488e-2# in zJ/K # 1.3806488e-23 # in J/K
hbar =6.582e-16 # in inverse sec #6.625e-34 # in J/s
z = ns * coeff
#r_1 = 0.5e-9
#r_2 = 0.5e-9
r_1 = 0.5
r_2 = 0.5
t = np.linspace(0.0,10000.0,100)
#ls = 8e-10#[1e-6 , 1e-5 , 1e-4]
#ls = linspace(0.1e-9, 10e-9, 100)
ls = linspace(0.5,2.5, 100)
#ls = [4.0]
#thetas = linspace(0.0,np.pi/2.0,100)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
#thetas = linspace(np.pi/5.0,np.pi/2.0,100)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
thetas = [np.pi/5.0,np.pi/4.0,np.pi/3.0]#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
eiz_m = 1.0

#height= np.zeros(shape=(len(ls),len(thetas)))
aiz = np.zeros(len(eiz_z))
def a(perp, par):
    return 2.0*(perp-1.0)/((perp+1.0)*(par-1.0))
aiz = a(eiz_x,eiz_z)

lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']

A = np.empty(len(ns))
g = np.empty(len(ns))
A0 = np.empty(len(ns))
sum_A = np.empty(len(ls))
Sum_A = np.empty(len(ls))
sum_A_thetas = np.empty(len(thetas))
g_t = np.empty(len(thetas))
A_t = np.empty(len(thetas))
G = np.empty(len(ls))
added = np.empty(len(ls))
added_thetas = np.empty(len(thetas))

pl.figure()
for p in range(len(thetas)):
	for k in range(len(ls)):
		for j in range(len(ns)):
			#A = ((eiz_z-1)/1)*((eiz_z-1)/1)*(2*((1 +3 * aiz)*(1+3 * aiz))+(1 - aiz)(1 - aiz)*cos(thetas[p]))
			#sum_A[k] = sum(A)
			#sum_A_thetas[p] = sum(A[j])

		#	A = np.nan_to_num(A)
			#A = [(3.0/8.0)*((eiz_zs-1.0)/1.0)*((eiz_zs-1.0)/1.0)*(2.0*((1.+3.*aizs)*(1.+3.*aizs))+(1.-aizs)*(1.0 - aizs)*cos(thetas[p])) for eiz_zs,aizs in zip(eiz_z,aiz)]
			g[j] = 2.0*((1.0+3.0*aiz[j])*(1.0+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p])
			A[j] = (3.0/8.0)*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*(2.0*((1.+3.0*aiz[j])*(1.+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p]))
		g_t[p] = sum(g[j])#+= 2.0*((1.0+3.0*aiz[j])*(1.0+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p])
		A_t[p] = sum(A[j])#+= 2.0*((1.0+3.0*aiz[j])*(1.0+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p])
		Sum_A[k] +=    (3.0/8.0)*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*(2.0*((1.0+3.0*aiz[j])*(1.0+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p]))
		Sum_A[k] = np.nan_to_num(Sum_A[k])
		sum_A[k] = sum(A[j])
		#pl.plot(ls,Sum_A, color = colors[p])#, linestyle = lts[k])	
			#print sum_A
			#sum_A_thetas[p] = sum(A[j])
		#G = -(1.6e-19*(kb*T)*(np.pi*np.pi*r_1*r_1*r_2*r_2)/(64.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A
		#G[k] = -((kb_J*T)*(np.pi*np.pi*r_1*r_1*r_2*r_2)/(64.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A
		G[k] = -((kb_J*T)*(np.pi*np.pi*r_1*r_1*r_2*r_2)/(64.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*A_t[p]
	#pl.plot(ns,g, color = colors[p])#, linestyle = lts[k])	
	#pl.plot(ls,G, color = colors[p])#, linestyle = lts[k])	
	#pl.semilogy(ls,-G, color = colors[p])#, linestyle = lts[k])	
#	pl.figure()

	pl.plot(ls,-G, color = colors[p])#, linestyle = lts[k])	
pl.xlabel(r'$l\,\,\,   [nm]$', size = 22)
pl.ylabel(r'$ - G(l,\theta)\,\,[zJ]$', size = 22)
pl.title(r'G(l, $\theta$) Non retarded skewed cylinders', size = 16)
pl.savefig('plots/131027_Hopkins_nonret_G_l_theta.pdf')

#pl.figure()
#pl.plot(ls,A_t/A_t[0])#, color = colors[p])#, linestyle = lts[k])	
#pl.xlabel(r'$\theta\,\,\,[radians]$', size = 22)
#pl.ylabel(r'$ A(\theta)/A(0)$', size = 22)
#pl.title(r'$\Sigma \Delta_{1} \Delta_{2} \times \rm{g(l = 4 nm, \theta)}$ for non retarded skewed cylinders', size = 16)
#pl.savefig('plots/131027_Hopkins_nonret_sum_theta.pdf')


	#pl.text(0.1,-0.15e-8,r'$\theta = %.2f \, radians$' % thetas[0], color='r')
	#pl.text(0.1,-0.30e-8,r'$\theta = %.2f \, radians$' % thetas[1], color='g')
	#pl.text(0.1,-0.45e-8,r'$\theta = %.2f \, radians$' % thetas[2], color='b')
#pl.legend()
#pl.savefig('plots/131010_Hopkins_CG10_eps2_x_z.png', dpi = 300 )
pl.show()
pl.close()






