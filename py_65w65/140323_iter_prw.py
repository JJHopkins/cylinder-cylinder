#!/usr/bin/env python
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as pl

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') # len(n), LDS, ~material response in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # len(n), LDS, ~material response in parallel direction
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # len(n), LDS, ~response of water, intervening medium

# Constants
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
Temp = 297              # [K] 
kbT = Temp * 1.3807e-23 # [J]

# Matsubara frequencies
ns = np.arange(0.,500.)         # index for Matsubara sum
zs = ns * coeff                 # thermal frequencies, they are multiples of n
Ls = np.arange(1e-9,1e-6,1e-9)  # separation distance between 2 cyclinders
#Ls = np.arange(1.e-9,1.e-6,10e-9) #np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders
t = np.linspace(0.,1e5,1e4)
u = np.linspace(0.,1e5,1e4) 
#u = np.linspace(0.,1e4,1e4) 
y = np.linspace(1.000001,100001.000001,100)

# Define functions
def Aiz(perp, par,med):
	return (2.0*med*(perp-med)*(par-med))/((perp + med)*(par*par -2.*par*med + med*med))
	#return (2.0*(perp-med)*med)/((perp+med)*(par-med))
    #multiplies aiz by 1 = (par-med)/(par-med) to get zero instead of
    #divergence

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A = np.zeros(shape = (len(Ls),len(ns)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

T,Y = np.meshgrid(t,y)

U,Y0 = np.meshgrid(u,y)
		
# Calculate N=0 term        
f_term0 = (1./np.sqrt(Y0 * Y0 -1.0)) * U*U*U*U\
        *np.exp(-2.* Y0 *U)\
        *(2.*(1. + 3.*a[0])*(1. + 3.*a[0])\
        +(1.-a[0])*(1.-a[0]))
Ft_term0 = np.sum(f_term0 , axis = 1)
Fty_term0 = np.sum(Ft_term0)

# Calculate 1 < N < 500 terms
for i,L in enumerate(Ls):
    print L
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A
        f = (1./np.sqrt(Y*Y -1.0))\
                * np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1.))\
                * (T / np.sqrt(T*T + 1.))\
                * (2. * (1. + 3.*a[j]) * (1. + 3.*a[j]) * T*T*T*T\
                + 4. * (1. + 2.*a[j] + 2.*a[j] + 3.*a[j] * a[j]) *T*T\
                + 4.*(1.+a[j])*(1.+a[j])\
                + (T*T*T*T + 4.*T*T + 4.)*(1.-a[j])*(1.-a[j]))
       
        Ft = np.sum(f , axis = 1)
        Fty = np.sum(Ft)
        # Matsubara terms 
        A[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Fty
        A[i,0] = (1./2) * delta[0]*delta[0]*Fty_term0
    sum_A = (kbT/(12.* np.pi)) * np.sum(A, axis = 1)
np.savetxt('65w65_prw_sum_A.txt',sum_A)
np.savetxt('65w65_prw_Ls.txt',Ls)

pl.figure()
pl.loglog(Ls,sum_A,label=r'$\mathcal{A}(\ell= %1.1f nm)= %3.2f \,\,\, \mathcal{A}(\ell= %1.1f nm)= %3.2f $'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[3],1e21*sum_A[3]))
pl.legend(loc = 'best')
pl.show()
#
#
#########################################
#for i,L in enumerate(Ls):
#	for j,n in enumerate(ns):
#		p[i,j] = Pn(eiz_w[j],zs[j],L,c)
#		f = (1./np.sqrt(Y**2 -1.0))*np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1))*(T/np.sqrt(T*T+1))*(2.*(1.-3.*a[j])*(1.-3.*a[j])*T**4 + 4.*(1.+2.*a[j]+2.*a[j]+3*a[j]*a[j])*T**2 + 4.*(1.+a[j])*(1.+a[j]) +  (T**4 + 4.*T**2+4.)*(1.-a[j])*(1.-a[j]) )
#	
#		Ft = np.sum(f,axis = 1)
#		Fty = np.sum(Ft)#, axis = 1)
#	#	print j,Fty
#		A[i,j] = delta[j]*delta[j]*(p[i,j]**5)*Fty
#		print j,A
#		sum_A = (kbT/(12.*np.pi)) * np.sum(A, axis = 1)
#np.savetxt('140316_65w65_prw_sum_A.txt',sum_A)
#np.savetxt('140316_65w65_prw_Ls.txt',Ls)
#
#pl.figure()
#pl.loglog(Ls,sum_A)
#pl.show()
