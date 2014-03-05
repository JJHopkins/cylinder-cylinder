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
c = 3.0*1e8*1e9
#r_1 = 1.0e-9
#r_2 = 1.0e-9
r_1 = 0.5
r_2 = 0.5
t = np.linspace(0.0,100.0,100)
ls = [4.0]#linspace(1.0e-9, 100.0e-9 , 100)
thetas = linspace(np.pi/100.0,np.pi/2.0,50)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
eiz_m = 1.0

#height= np.zeros(shape=(len(ls),len(thetas)))
aiz = np.zeros(len(eiz_y))
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
sum_A_2 = np.empty(len(ls))
G = np.empty(len(thetas))
G_2 = np.empty(len(thetas))
added = np.empty(len(thetas))

pl.figure()
for p in range(len(thetas)):
	for k in range(len(ls)):
		for j in range(len(ns)):
			int_g[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(2.0*((1.+3.*aiz[j])*(1.+3.*aiz[j])*t*t*t*t+2.*(1.+2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])*t*t+2.0*(1.0+aiz[j])*(1.0 + aiz[j]))),t)

			int_g_2[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(((1.0 - aiz[j])*(1.0 - aiz[j])*(t*t + 2.0)*(t*t + 2.0))),t)

			A[j] = (kb)*T/32.0*((eiz_x[j]-1.0)/1.0)*((eiz_x[j]-1.0)/1.0)*\
			1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g[j]

			A_2[j] = (kb)*T/32.0*((eiz_x[j]-1.0)/1.0)*\
			((eiz_x[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g_2[j]	
			
			sum_A = sum(A[j])

			sum_A_2 = sum(A_2[j])

			G = -((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A

			G_2 = -((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_2
	added[p] = -1.0/(kb*T)*(-((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A)+\
	(-((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/\
	(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_2)	
	pl.plot(thetas,added)#, color = colors[p])#, linestyle = lts[k])	
#pl.axis([0.0,0.003,-0.1e-40,0.0])
pl.show()





