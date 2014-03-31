#!/usr/bin/env python

import numpy as np
from scipy.integrate import quad

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') # len(n), LDS, ~material response in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # len(n), LDS, ~material response in parallel direction
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # len(n), LDS, ~response of water, intervening medium

# Constants
c = 2.99e8   # light speed in m/s
coeff = 2.411e14 # in rad/s
Temp = 297 # K 
kbT = Temp * 1.3807e-23 # in J

# Matsubara frequencies
ns = np.arange(0.,500.)  # index for Matsubara sum
zs = ns * coeff           # thermal frequencies, they are multiples of n

Ls0 = 4.e-9 
Ls1 = 4.e-8 
t = np.linspace(0.,100000.,1000)
y = np.linspace(1.000001,100001.000001,1000)

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l,v):
	return np.sqrt(e)*zn*l*(1./v)

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)
p0 = Pn(eiz_w,zs,Ls0,c)
p1 = Pn(eiz_w,zs,Ls1,c)

T,Y = np.meshgrid(t,y)
#f = np.zeros(len(ns))
A0 = np.empty(len(ns))
A1 = np.empty(len(ns))

for j,n in enumerate(ns):
	f0 = (1./np.sqrt(Y**2 -1.0))*np.exp(-2.*Y*p0[j]*np.sqrt(T*T+1))*(T/np.sqrt(T*T+1))*(2.*(1.-3.*a[j])*(1.-3.*a[j])*T**4 + 4.*(1.+2.*a[j]+2.*a[j]+3*a[j]*a[j])*T**2 + 4.*(1.+a[j])*(1.+a[j]) +  (T**4 + 4.*T**2+4.)*(1.-a[j])*(1.-a[j]) )

	f1 = (1./np.sqrt(Y**2 -1.0))*np.exp(-2.*Y*p1[j]*np.sqrt(T*T+1))*(T/np.sqrt(T*T+1))*(2.*(1.-3.*a[j])*(1.-3.*a[j])*T**4 + 4.*(1.+2.*a[j]+2.*a[j]+3*a[j]*a[j])*T**2 + 4.*(1.+a[j])*(1.+a[j]) +  (T**4 + 4.*T**2+4.)*(1.-a[j])*(1.-a[j]) )

	Ft0 = np.sum(f0,axis = 1)
	Ft1 = np.sum(f1,axis = 1)

	Fty0 = np.sum(Ft0)
	Fty1 = np.sum(Ft1)

	print j,Fty0
	print Fty1

	A0[j] = delta[j]*delta[j]*(p0[j]**5)*Fty0
	A1[j] = delta[j]*delta[j]*(p1[j]**5)*Fty1
	#print A

	sum_A0 = (kbT/(12.*np.pi)) * np.sum(A0)
	sum_A1 = (kbT/(12.*np.pi)) * np.sum(A1)

print sum_A0
print sum_A1



