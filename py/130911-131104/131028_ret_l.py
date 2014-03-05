#!/usr/bin/python

import numpy as np
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import axis as ax
import matplotlib               
from matplotlib import pyplot as pl
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

x_t_eV ,y_t  = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x_eV ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y_eV ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z_eV ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

x_t =x_t_eV# /6.582e-16
x_x =x_x_eV# /6.582e-16
x_y =x_y_eV# /6.582e-16
x_z =x_z_eV# /6.582e-16

eiz_x = np.loadtxt('data/eiz_x_output.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') # parallel,axial

coeff = 0.159 # in eV #(2.41*1e14) # in rad/s
#coeff = 2.41 # in rad/s
ns = np.arange(1.0,501)
T = 300.0 # r.t. in K
#kb_J =1.3806488e-23# in zJ/K # 1.3806488e-23 # in J/K
kb_J =1.3806488# in zJ/K # 1.3806488e-23 # in J/K
z = ns * coeff
r_1 = 0.5e-9
r_2 = 0.5e-9
c = 2.99e8

#ls = np.linspace(1e-9, 1000e-9, 10)
ls = np.linspace(1e-9, 1000e-9, 10)
thetas = np.array([np.pi/5.0])
eiz_m = 1.0

lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']

aiz = []
sum_A = np.empty(len(ls))
sum_A_2 = np.empty(len(ls))

G = np.empty(len(thetas))

def Aiz(perp, par):
	return 2.0*(perp-1.0)/((perp+1.0)*(par-1.0))

def ys(a,time,L, N):
	term0 = time    / (time*time+1.0) 
	term1 =	time**4 * 2.0*(1. + 3.*a)*(1.+3.*a)
	term2 = time**2 * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a)
	term3 =           4.0*(1. + a)*(1.0 + a) 
	term4 = np.exp(1e8*(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8))
	#term4 =	1e10*(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8) +\
	#	1e10*(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**2/2.0 +\
	#	1e10*(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**3/6.0
	#term5 = 1.0 + 1e-10*term4
	#print 'ys term0', term0
	#print 'ys term1', term1
	#print 'ys term2', term2
	#print 'ys term3', term3
	#print 'ys term4', term4
	#print '----'
	return term0 * term1 * term2 * term3 * term4

def y_2s(a,time, L, N): 
	term0 = time    / (time*time+1.0) 
	term1 = (1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0)                     
	term2 = np.exp(1e8*(-2.0*L*coeff*N*np.sqrt(time*time + 1.0)/2.99e8))
	#term2 =	1e10*(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8) +\
	#	1e10*(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**2/2.0 +\
	#	1e10*(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**3/6.0
	#term2 = 1e12*(-2.0*L*coeff*N*np.sqrt(time*time + 1.0)/2.99e8)
	#print 'ys term0', term0
	#print 'ys term1', term1
	#print 'ys term2', term2
	#print '----'
	return term0 * term1 * term2


def As(eizz,L,N,Y): 
	term0 = (1.0/32)*kb_J*T
	term1 = ((eizz-1.0)/1.0)*((eizz-1.0)/1.0)
	term2 = 1.0 * (coeff*N)**4 * L**4 * (2.99)**(-4)
	term3 = Y
	
	return term0 * term1 * term2 * term3

def A_2s(eizz, L , N ,Y):		
	term0 = (1.0/32)*kb_J*T
	term1 = ((eizz-1.0)/1.0)*((eizz-1.0)/1.0)
	term2 = 1.0 * (coeff*N)**4 * L**4 * (2.99)**(-4)
	term3 = Y
	
	return term0 * term1 * term2 * term3

aiz = Aiz(eiz_z,eiz_x) # of length = len(ns)

dt = 1
#ts = np.arange(1e-9,101e-9,dt)
ts = np.arange(1,101,dt)

A   = np.zeros(shape=(len(ns),len(ls)))
A_2 = np.zeros(shape=(len(ns),len(ls)))

pl.figure()
from pylab import pause
for j,n in enumerate(ns):

	for k,l in enumerate(ls):

		# Integrand:
		y_arg   = ys(aiz[j],ts,l,n)
		y_2_arg = y_2s(aiz[j],ts,l,n)
		

		# Integral:
		y   = trapz(y_arg,ts,dt)
		y_2 = trapz(y_2_arg,ts,dt)

		#print As(eiz_z[j],l,n,y)
		#print A_2s(eiz_z[j],l,n,y_2)

		A[j,k]   = As(eiz_z[j],l,n,y)
		A_2[j,k] = A_2s(eiz_z[j],l,n,y_2) 		


#pl.show()
print '----'

EL = ls**(-4)
print 'EL =', EL
sum_A = np.sum(A,axis=0)
print 'sum_A = ', sum_A
sum_A_2 = np.sum(A_2,axis=0)
print 'sum_A_2 = ', sum_A_2
G = -(EL*np.pi*r_1*r_1*r_2*r_2)*(sum_A + sum_A_2*np.cos(2.0*np.pi/5))/(2.0*np.sin(np.pi/5)) 
print G
#G=np.zeros(shape=(len(ls),len(thetas)))
#for t,theta in enumerate(thetas):
    #G[:,t] = -((L*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.sin(theta)))\
#	*(sum_A + sum_A_2*np.cos(2.0*theta)) 
#G = -1.0*L*(sum_A + sum_A_2 )

pl.figure()
pl.loglog(ls,-G)#,color = colors[j])
pl.show()
