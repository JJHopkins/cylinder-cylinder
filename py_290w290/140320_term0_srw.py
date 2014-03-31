#!/usr/bin/env python
import numpy as np
from scipy.integrate import trapz
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
T = np.linspace(0.,1e4,1e4)
U = np.linspace(0.,1e4,1e4) 

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
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)
for i,L in enumerate(Ls):
    print L
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))*(T/(T*T+1.))\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))
        # A0(n=0) term 
        f0_term0 = U*U*U * np.exp(-2.* U)\
                *2.*(1.+3.*a[0])*(1.+3.*a[0])
        # Integrand A2                
        f2 = np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                *(T/(T*T+1.))\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))
        # A2(n=0) term
        f2_term0 =  U*U*U * np.exp(-2.* U)\
                *(1.-a[0])*(1.-a[0])

        Ft0 = np.sum(f0)
        Ft0_term0 = np.sum(f0_term0)
        Ft2 = np.sum(f2)
        Ft2_term0 = np.sum(f2_term0)
        # Matsubara terms 
        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
np.savetxt('290w290_srw_sum_A0.txt',sum_A0)
np.savetxt('290w290_srw_sum_A2.txt',sum_A2)
np.savetxt('290w290_srw_Ls.txt',Ls)

pl.figure()
pl.loglog(Ls,sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f\,\,\,\        mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
pl.loglog(Ls,sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=%1.1f nm)=%3.2f\,\,\,\mathcal{A^{(2)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[3],1e21*sum_A2[3]))
pl.legend(loc = 'best')
pl.show()
