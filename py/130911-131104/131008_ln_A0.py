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
t = np.linspace(1.0,1e17,100)
#ls = 8e-10#[1e-6 , 1e-5 , 1e-4]
ls = linspace(0.1 , 4.0 , 50)
thetas = linspace(np.pi/4.0,(7./8)*np.pi,6)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
eiz_m = 1.0

#height= np.zeros(shape=(len(ls),len(thetas)))
aiz = np.zeros(len(eiz_z))
def a(perp, par):
    return 2.0*(perp-1.0)/((perp+1.0)*(par-1.0))
aiz = a(eiz_x,eiz_z)

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
pre_int_g = np.empty(len(ls))

pl.figure()
for p in range(len(thetas)):
	for k in range(len(ls)):
		for j in range(len(ns)):

			#pre_int_g[j] = np.log(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			#(2.0*((1.+3.*aiz[j])*(1.+3.*aiz[j])*t*t*t*t+2.*(1.+2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])*t*t+2.0*(1.0+aiz[j])*(1.0 + aiz[j]))))
			#int_g[j] = np.exp(pre_int_g)
			int_g[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(2.0*((1.+3.*aiz[j])*(1.+3.*aiz[j])*t*t*t*t+2.*(1.+2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])*t*t+2.0*(1.0+aiz[j])*(1.0 + aiz[j]))),t)

			A[j] = 1.0/32.0*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*\
			1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g[j]

			sum_A[k] = sum(A[j])
		print sum_A	
	G = -((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A

	savetxt("data/G0_output.txt", G)
	pl.plot(ls,G, color = colors[p])#, linestyle = lts[k])	
pl.show()




