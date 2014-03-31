#!/usr/bin/env python
#
from __future__ import division
import sys 
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as pl
from matplotlib import axis as ax
from pylab import show
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from pylab import pause
from scipy.integrate import dblquad
from matplotlib.backends.backend_pdf import PdfPages
from scipy import integrate

# We want to solve for the fully retarded interaction free energy, G, between two identical, 
# anisotropic, parallel cylinders in water.
# The inputs are the dielectric responses of the cylinder material and water as functions of n
# For each separation (length) and value of n, we want to find the double integral over dh and dt.
# We then sum over the values of n, this gives us a hamaker coefficent (A) for each value of 
# separation (length).

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') # len(n), LDS, ~material response in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # len(n), LDS, ~material response in parallel direction
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # len(n), LDS, ~response of water, intervening medium

# Constants
r_1 = 1.0e-9 # radius of cylinder_1
r_2 = 1.0e-9 # radius of cylinder_2
c = 2.99e8   # light speed in m/s
#c = 1.0 #2.99e8   # light speed in m/s
coeff = 2.411e14 # in rad/s
#coeff = 0.159

# Matsubara frequencies
ns = np.arange(0.,500.)  # index for Matsubara sum
z = ns * coeff           # thermal frequencies, they are multiples of n
ls = np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders

def Aiz(perp, par,med):
	'''
	response mismatch ratio, uses differences between cylinder and water responses (eiz)
	'''
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def integrand(t, N, h, a, eizw, L):
    	"""
    	function of time, height, aiz(n), eiz_w(n), length, and n. Needs to be integrated over
	dtdy for each value of length and for aiz,eiz,and n. 
   	"""
	term0 = 1./ np.sqrt(h*h - 1.0 )
	term1 = ( t    / np.sqrt(t*t+1.0)     )
	#term2 = (-2.0 * h * np.sqrt(eizw)* L * coeff * N * (1./c) * np.sqrt(t*t + 1.0))
	term2 = (-2.0 * h * np.sqrt(eizw)* L * coeff * N* (1./(6.58e-16)) * (1./c) * np.sqrt(t*t + 1.0))
	term3 =	( t**4 * 2.0*(1. + 3.*a)*(1.+3.*a)     )
	term4 = ( t**2 * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a))
	term5 = (          4.0*(1. + a)*(1.0 + a)         )
	term6 = ((1.- a)*(1.- a)*(t * t  + 2.0)*(t * t + 2.0))                     
	print 'g0 term0', term0
	print 'g0 term1', term1
	print 'g0 term2', term2
	print 'g0 exp(term2)', np.exp(term2)
	#print 'g0 exp(term2)^pwr', (np.exp(term2))**(coeff/c)
	print 'g0 term3', term3
        print 'g0 term4', term4
        print 'g0 term5', term5
	print 'g0 term6', term6
	print '----'
	return term0 * term1 * (np.exp(term2)**(coeff/c))*( term3 + term4 + term5 + term6)

def inner_integral(N, h , a, eizw, L):
	r"""
	first integral is over dt from 0 to inf.  
	"""
	return integrate.quad(integrand, 3.3e-17, 0.9, args=(N,h, a,eizw,L))[0]

def As(eizz,eizw,L,N,Integral): 
	'''
	Hamaker coefficent, gives strength of vdW interaction
	'''
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Integral * eizw **(5./2) * (coeff*N)**5 * L**5 / (c**5))
	return  term1 *  term2

# Calculate anisotropy metric using dielectric mismatches
aiz = Aiz(eiz_x,eiz_z, eiz_w) 

# Intialize
#I   = np.zeros(shape=(len(ls),len(ns)))  # double integral result 
A   = np.zeros(shape=(len(ls),len(ns)))  # hamaker coefficient
EL = np.zeros(len(ls))                   # prefactor: length^(-5)
G_l_t_dt = np.zeros(len(ls))             # interaction free energy per unit length   


for j,length in enumerate(ls):

	sum_A = np.empty(len(ls))
	sum_A_2 = np.empty(len(ls))

	for k,n in enumerate(ns):
	
		I = integrate.quad(lambda h: inner_integral(n,h,aiz[k],eiz_w[k],length), 1.1, 1.0e9)#np.inf)
		print I[0] 
		
		# plug in the intergral result for each length,n pair

		A[j,k]   =   As(eiz_z[k],eiz_w[k],length,n,I[0])
	A[j,0]   *= 0.5  
	sum_A = np.sum(A,axis=1)

for k,length in enumerate(ls):
	EL[j] = 1./(length*length*length*length*length)
	G_l_t_dt[j] = (1.602e-19 / 4.11e-21) * (1./32) * EL[j]*r_1*r_1*r_2*r_2*(1./(12*np.pi))*(sum_A[j])

