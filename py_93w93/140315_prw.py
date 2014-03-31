#/usr/bin/env python

import numpy as np
#import scipy.integrate as integrate
from scipy.integrate import quad

#def integrand(t,n,x,b):
#    return (1./np.sqrt(x*x - 1.)) * t * (1./np.sqrt(t*t + 1.)) * np.exp(-2. * x * np.sqrt(t*t +1.)) 
#
#def expint(n,x,b):
#    return quad(integrand, 1, np.inf, args=(n, x,b))[0]
#
#result = quad(lambda x: expint(n,x,b), 0, np.inf)

#ns = np.arange(1.,5.) 
#bs = np.arange(11.,15.) 
#
#for j,n in enumerate(ns):
#	for k,b in enumerate(bs):
#		result = quad(lambda x: expint(n, x,b), 0, np.inf)

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') # len(n), LDS, ~material response in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # len(n), LDS, ~material response in parallel direction
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # len(n), LDS, ~response of water, intervening medium

# Constants
c = 2.99e8   # light speed in m/s
#c = 1.0 #2.99e8   # light speed in m/s
coeff = 2.411e14 # in rad/s
#coeff = 0.159
Temp = 297 # K 
kbT = Temp * 1.3807e-23 # in J
#kbT = Temp * 8.6173e-5  # in eV
# Matsubara frequencies
ns = np.arange(0.,500.)  # index for Matsubara sum
zs = ns * coeff           # thermal frequencies, they are multiples of n
Ls = [4.e-9, 4.e-8, 4.e-7, 4.e-6] #np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders
#Ls = 4.e-9 #np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders

# 1 eV = 1.602e-19 J = 0.016 zJ
eVtoJ =  1.602e-23
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

#def Pn(e,zn,l,v):
#	return np.sqrt(e)*zn*l*(1./v)

def Pn_L(e,zn):
	return np.sqrt(e)*zn*(1./c)

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)
pn_L = Pn_L(eiz_w,zs)

#T = np.linspace(0.,100000.,5)
#Y = np.linspace(1.000001,100001.000001,6)
#T,Y = np.meshgrid(t,y)
#f = np.zeros(shape = (len(ns),len(Ls),len(Ls),len(Ls)))
A = np.zeros(shape = (len(Ls),len(ns)))

for i,L in (enumerate(Ls)):
	for j,n in enumerate(ns):
		T = np.linspace(0.,100000.,5)
		Y = np.linspace(1.000001,100001.000001,5)
		f = (1./np.sqrt(Y**2 -1.0))*np.exp(-2.*Y*pn_L[j]*L*np.sqrt(T*T+1)) \
		* (T/np.sqrt(T*T+1))*(2.*(1.-3.*a[j])*(1.-3.*a[j])*T**4 \
		+ 4.*(1.+2.*a[j]+2.*a[j]+3*a[j]*a[j])*T**2 + 4.*(1.+a[j])*(1.+a[j]) \
		+  (T**4 + 4.*T**2+4.)*(1.-a[j])*(1.-a[j]) )# for L in Ls]
		
		Ft = np.sum(f,axis = 1)
		Fty = np.sum(Ft)#, axis = 1)
		#print j,Fty
		#A[j] = delta[j]*delta[j]*(p[j]**5)*Fty
		
		A = delta[j]*delta[j]*(pn_L[j]**5)*(L**5)*Fty#[i]# for F in Fty]
		sum_A[i] = (kbT/(12.*np.pi)) * np.sum(A, axis = 1)
		
		#p = [Pn(eiz_w,zs,L,c) for L in Ls]
		#f = [(1./np.sqrt(Y**2 -1.0))*np.exp(-2.*Y*p[j]*np.sqrt(T*T+1)) \
		#* (T/np.sqrt(T*T+1))*(2.*(1.-3.*a[j])*(1.-3.*a[j])*T**4 \
		#+ 4.*(1.+2.*a[j]+2.*a[j]+3*a[j]*a[j])*T**2 + 4.*(1.+a[j])*(1.+a[j]) \
		#+  (T**4 + 4.*T**2+4.)*(1.-a[j])*(1.-a[j]) )
	
#print A
#A[0] = 0.5*A[1]
#sum_A = (kbT/(12.*np.pi)) * np.sum(A, axis = 1)
#sum_A = (kbT/(12.*np.pi)) * np.sum(A)
#print sum_A
#print sum_A
# h_bar_eV = 6.5821e-16 eVs
# h_bar = 1. #1.0546e-34 #in Js


