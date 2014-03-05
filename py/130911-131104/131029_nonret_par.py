#!/usr/bin/python

import numpy as np
from numpy import log,exp
#from scipy import *
import scipy.optimize as opt
from scipy.integrate import trapz
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
x_w_eV ,y_w = np.loadtxt('data/LO_water.csv',delimiter=',',unpack=True, usecols = [0,1])
x_t =x_t_eV /6.582e-16
x_x =x_x_eV /6.582e-16
x_y =x_y_eV /6.582e-16
x_z =x_z_eV /6.582e-16
x_w =x_w_eV /6.582e-16
eiz_x = np.loadtxt('data/eiz_x_output.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') # parallel,axial
eiz_w = np.loadtxt('data/eiz_w_output.txt') # parallel,axial

coeff = 2.41 # in rad/s
ns = arange(1.0,501)
T = 300.0 # r.t. in K
kb_J =1.3806488e-23# in zJ/K # 1.3806488e-23 # in J/K
hbar =6.582e-16 # in inverse sec #6.625e-34 # in J/s
z = ns * coeff
r_1 = 0.5e-9
r_2 = 0.5e-9
#r_1 = 0.5
#r_2 = 0.5
c = 2.99e8

#ls = linspace(0.1e-9, 10e-9, 100)
ls = linspace(1.1e-9, 1000e-9, 10)
#ls = [4.0]
#thetas = linspace(0.0,np.pi/2.0,100)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
#thetas = [np.pi/5.0,np.pi/4.0,np.pi/3.0]#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
thetas = [0.0]#np.pi/5.0,np.pi/4.0,np.pi/3.0]#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#

#height= np.zeros(shape=(len(ls),len(thetas)))
lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']
#
aiz = []#np.zeros(len(eiz_z))
#A = np.empty(len(ns))
#A_2 = np.empty(len(ns))
#y = np.empty(len(ns))
#y_2 = np.empty(len(ns))
#g = np.empty(len(ns))
#g_int = np.empty(len(ns))
#A0 = np.empty(len(ns))
#w = np.empty(len(ls))
sum_A = np.empty(len(ls))
sum_A_2 = np.empty(len(ls))
#sum_A_thetas = np.empty(len(thetas))
#g_t = np.empty(len(thetas))
#A_t = np.empty(len(thetas))
#G = np.empty(len(ls))
#added = np.empty(len(ls))
#added_thetas = np.empty(len(thetas))
#
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def ys(a): 
	return (3. + 5.*a + 5.*a + 19.*a*a)
#def ys(a,time,L, med,N):
#	return (\
#		(((time**5/(time*time+1.0)) * 2.0*(1. + 3.*a)*(1.+3.*a)     ) + \
#		((time**3/(time*time+1.0)) * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a)) + \
#		((time   /(time*time+1.0)) * 4.0*(1. + a)*(1.0 + a))) * \
#		np.exp((-2.0)*L*med**(1./2)*(1./(2.99e8))*coeff*N*(time*time + 1.0)**(1./2)))
#def y_2s(a,time, L,med, N): 
#	return (\
#		((time / (time * time + 1.0)) * (1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0)) * \
#		np.exp((-2.0)*L*med**(1./2)*(1./(2.99e8))*coeff*N*(time*time  + 1.0)**(1./2)))
def As(eizz, eizw,Y): 
	return (((eizz-eizw)/eizw)*((eizz-eizw)/eizw)*Y)
		#np.exp(Y)
aiz = np.zeros(len(eiz_z))

for i,ezz in enumerate(eiz_z):
	aiz = [Aiz(ezx,ezz,ezw) for ezx,ezw in zip(eiz_x,eiz_w)]
	print aiz[i],eiz_w[i],i

print '****'
ts = linspace(1e-9,1000e-1)
dt =10e-9
L = np.zeros(len(ls))
#for j,n in enumerate(ns):
for p, l in enumerate(ls):
	L[p] = l**(-5)
	#print thetas,l
	y   = [ys(az) for az in aiz]
	A   = [As(ezz,ezw,yy) for ezz,ezw,yy in zip(eiz_z,eiz_w,y)]
	print A[499]
	sum_A[p]= sum(A)
	G = -(L*(9.0*pi*r_1*r_1*r_2*r_2)/(64.0*32.))*(sum_A) 
pl.figure()
pl.loglog(ls,-G)#,color = colors[j])
pl.show()
