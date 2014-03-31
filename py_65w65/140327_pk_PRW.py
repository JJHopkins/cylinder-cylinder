#!/usr/bin/env python
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # LDS in parallel direction
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # LDS of water, intervening medium
eiz_w[0] = 79.0

# Constants
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
Temp = 297              # [K] 
kbT = Temp * 1.3807e-23 # [J]

# Matsubara frequencies
ns = np.arange(0.,500.) 
zs = ns * coeff         

# Intersurface eparation distance between 2 cyclinders
Ls = np.arange(1e-9,1e-7,10e-9)        

#Integration vars
t  = np.linspace(0.,2.**9, 1.+2.**9)
u  = np.linspace(0.,2.**9, 1.+2.**9)
y  = np.linspace((1.+1e-5),2.**9, 1.+2.**9)
y0 = np.linspace((1.+1e-5),2.**9, 1.+2.**9)

# Vectorize
T,Y = np.meshgrid(t,y)
U,Y0 = np.meshgrid(u,y0)
		
# Initialize
p = np.zeros(shape = (len(Ls),len(ns)))
A = np.zeros(shape = (len(Ls),len(ns)))

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

# Calculate N=0 term        
# Integrand
f_term0 = (1./(Y0*np.sqrt(1.0 - (1./Y0)))) *U*U*U*U\
        *np.exp(-2.* Y0 *U)\
        *(2.*(1. + 3.*a[0])*(1. + 3.*a[0])\
        +(1.-a[0])*(1.-a[0]))
# Double Integral
Ft_term0 = romb(f_term0, axis = 1)
Fty_term0 = romb(Ft_term0, axis = 0)

# Calculate 1 < N < 500 terms
for i,L in enumerate(Ls):
    print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand
        f = (1./(Y*np.sqrt(1.0 - (1./Y))))\
                * np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1.))\
                * (T / np.sqrt(T*T + 1.))\
                * (2. * (1. + 3.*a[j]) * (1. + 3.*a[j]) * T*T*T*T\
                + 4. * (1. + 2.*a[j] + 2.*a[j] + 3.*a[j] * a[j]) *T*T\
                + 4.*(1.+a[j])*(1.+a[j])\
                + (T*T*T*T + 4.*T*T + 4.)*(1.-a[j])*(1.-a[j]))
        # Double integral
        Ft = romb(f , axis = 1)
        Fty =romb(Ft, axis = 0)
    A[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Fty
    A[i,0] = (1./2) * delta[0]*delta[0]*Fty_term0
sum_A = (kbT/(12.* np.pi)) * np.sum(A, axis = 1)

print 'A(separation) = ',sum_A
print 'Contribution to A from n=0 term = ', (kbT/(12.*np.pi))*A[:,0]

np.savetxt('data/65w65_npk_A_parallel_ret.txt',sum_A)
np.savetxt('data/65w65_npk_Lengths_parallel_ret.txt',Ls)

# PLOTS
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

# Log-log
pl.figure()
pl.loglog(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_par)
pl.title(title('6','5','parallel'))
pl.legend(loc = 'best')
pl.savefig(svfig('65pk','65','loglog_parallel'))
pl.show()

# Semilog
pl.figure()
pl.semilogy(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.title(title('6','5','parallel'))
pl.legend(loc = 'best')
pl.axis([0.0,500, 1e1,1e3])
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig(svfig('65pk','65','parallel'))
pl.legend(loc = 'best')
pl.show()
