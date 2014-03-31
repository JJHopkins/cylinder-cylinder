#!/usr/bin/env python
import numpy as np
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

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])
Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    print L
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))*(T/(T*T+1.))\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))
        # Integrand A2                
        f2 = np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                *(T/(T*T+1.))\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))
        Ft0 = np.sum(f0)
        Ft2 = np.sum(f2)
        # Matsubara terms 
        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
np.savetxt('91w91_srw_sum_A0.txt',sum_A0)
np.savetxt('91w91_srw_sum_A2.txt',sum_A2)
np.savetxt('91w91_srw_Ls.txt',Ls)

A_py_par  = r'$\mathcal{A}_{\parallel}\sf{[python]}$'
A0_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
A2_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
A_GH_par  = r'$\mathcal{A}_{\parallel}\sf{[ G.H. ]}$'

x_ax = r'$\,\ell\,\,\,\rm{[nm]}$'
y_ax_par = r'$\mathrm{\mathcal{A}_\parallel(\ell)}\,\,\,\rm{[zJ]}$'
y_ax_per = r'$\mathrm{\mathcal{A}_\perp (\ell)}\,\,\,\rm{[zJ]}$'

def title(cnt1,cnt2,orientation):
	return r'$\mathrm{[%s,%s]\,\,Hamaker\,coeff:\, %s \,in\,water,\,retarded}$'%(cnt1,cnt2,orientation)

def svfig(cnt1,cnt2,orientation):
	return 'plots/140322_%sw%s_HCs_%s_ret.pdf'%(cnt1,cnt2,orientation)

pl.figure()
pl.loglog(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[3],1e21*sum_A0[3]))
pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[3],1e21*sum_A2[3]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.title(title('9','1','perpendicular'))
pl.legend(loc = 'best')
pl.savefig(svfig('9','1','loglog, perpendicular'))
pl.show()

pl.figure()
pl.semilogy(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[3],1e21*sum_A0[3]))
pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[3],1e21*sum_A2[3]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.title(title('9','1','perpendicular'))
pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig(svfig('9','1','perpendicular'))
pl.show()
pl.legend(loc = 'best')
pl.show()

