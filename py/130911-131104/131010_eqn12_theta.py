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

x_scl ,y_scal = np.loadtxt('data/Y26004L.csv',delimiter=',',unpack=True, usecols = [0,1])
x_t ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

eiz_x = np.loadtxt('data/eiz_x_output.txt') #parallel
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') #perpendicular

# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159e15 # in eV #(2.41*1e14) # in rad/s
#coeff = (2.41*1e14) # in rad/s
ns = arange(0.0,100.0)
z = ns * coeff
#kb =6.626e-34 # 8.6173e-5 #in eV/K
#kb = 8.6173e-5 #in eV/K
T = 300.0 # r.t. in K
c = 3.0e17
#r_1 = 1.0e-9
#r_2 = 1.0e-9
r_1 = 0.5
r_2 = 0.5
t = np.linspace(0.0,1000.0,100)
#ls = 8e-10#[1e-6 , 1e-5 , 1e-4]
#ls = linspace(0.1, 12.0, 100)
ls = [4.0]
thetas = linspace(np.pi/100.0,np.pi/2.0,50)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
#thetas = linspace(np.pi/3.0,np.pi/2.0,6)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
eiz_m = 1.0

#height= np.zeros(shape=(len(ls),len(thetas)))
aiz = np.zeros(len(eiz_z))
def a(perp, par):
    return 2.0*(perp-1.0)/((perp+1.0)*(par-1.0))
aiz = a(eiz_z,eiz_x)

lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']

int_g_arg = np.empty(len(ns))
int_g = np.empty(len(ns))
int_g_2 = np.empty(len(ns))
A = np.empty(len(ns))
A_2 = np.empty(len(ns))
A0 = np.empty(len(ns))
A2 = np.empty(len(ns))
sum_A = np.empty(len(ls))
sum_A_thetas = np.empty(len(thetas))
sum_A_2 = np.empty(len(ls))
sum_A_2_thetas = np.empty(len(thetas))
G = np.empty(len(ls))
G_2 = np.empty(len(ls))
added = np.empty(len(ls))
added_thetas = np.empty(len(thetas))

pl.figure()
for p in range(len(thetas)):
	for k in range(len(ls)):
		for j in range(len(ns)):
			int_g[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(2.0*((1.+3.*aiz[j])*(1.+3.*aiz[j])*t*t*t*t+2.*(1.+2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])*t*t+2.0*(1.0+aiz[j])*(1.0 + aiz[j]))),t)

			y = np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)

			int_g_2[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(((1.0 - aiz[j])*(1.0 - aiz[j])*(t*t + 2.0)*(t*t + 2.0))),t)

		#	A[j] = (kb)*T/32.0*((eiz_x[j]-1.0)/1.0)*((eiz_x[j]-1.0)/1.0)*\
		#	1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
		#	int_g[j]

		#	A_2[j] = (kb)*T/32.0*((eiz_x[j]-1.0)/1.0)*\
		#	((eiz_x[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
		#	int_g_2[j]	
			
			A[j] = 1.0/32.0*((eiz_x[j]-1.0)/1.0)*((eiz_x[j]-1.0)/1.0)*\
			1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g[j]

			A_2[j] = 1.0/32.0*((eiz_x[j]-1.0)/1.0)*\
			((eiz_x[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g_2[j]	
			
			sum_A[k] = sum(A[j])
			sum_A_thetas[p] = sum(A[j])
			
			sum_A_2[k] = sum(A_2[j])
			sum_A_2_thetas[p] = sum(A_2[j])

		added_thetas[p] = -(-((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_thetas[p])+\
		(-((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/\
		(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_2_thetas[p])	

		G = -((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A

		G_2 = -((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_2

		added = (-((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A)+(-((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/\
		(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_2)	

	#pl.plot(ls,added, color = colors[p])#, linestyle = lts[k])	


pl.figure()
pl.plot(thetas,added_thetas*10**6,label=r'$G(l=4.0\,nm,\theta)=-\frac{(\pi R_1^{2})(\pi R_2^{2})}{2\pi l^{4}\sin{\theta}}\left({\cal A}^{(0)}(l)+{\cal A}^{(2)}(l)\cos 2\theta\right)$') 
pl.xlabel(r'$\theta\,\,\,[radians]$', size = 24)
pl.ylabel(r'$G(l,\theta)\,\,\times\,10^{-6}$', size = 24)
pl.legend()
pl.title(r'Eqn 12 $G(l,\theta)\,with\,\varepsilon\prime\prime_{m}=1\,and \, l=4.0\,nm$')
#pl.savefig('plots/131010_Hopkins_CG10_eps2_x_z.png', dpi = 300 )
pl.savefig('plots/131010_Hopkins_egn12_G_vs_theta.pdf')
	#pl.axis([0.0,0.003,-0.1e-40,0.0])





