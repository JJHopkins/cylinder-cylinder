#!/usr/bin/env python
import numpy as np
from scipy.integrate import simps #quad
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
ns = np.arange(1.,501.)         # index for Matsubara sum
#Ns = np.arange(1.,500.)         # index for Matsubara sum
zs = ns * coeff                 # thermal frequencies, they are multiples of n

#Ls = np.arange(0.,1e2)#1e-9,101e-9,10e-9)  # separation distance between 2 cyclinders
Ls = np.linspace(1e-9,150e-9,15)  # separation distance between 2 cyclinders
#Ls = np.arange(1.e-9,1.e-6,10e-9) #np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders
t = np.linspace(0.,101,10003)#,100)
#u = np.linspace((1.+1e-5),1e2,99) 
u = np.linspace(0.,201,403)#,99) 
y = np.linspace((1.00001),102.00001,10003)#,98)
y0 = np.linspace((1.+1e-5),202,403)#,98)
#y = np.linspace(1.000001,1001.000001)

#Ls = np.arange(1e-9,1e-6,1e-9)  # separation distance between 2 cyclinders
##Ls = np.arange(1.e-9,1.e-6,10e-9) #np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders
#t = np.linspace(0.,1e5,1e4)
#u = np.linspace(0.,1e5,1e4) 
##u = np.linspace(0.,1e4,1e4) 
#y = np.linspace(1.000001,100001.000001,100)

# Define functions
def Aiz(perp, par,med):
	#return (2.0*med*(perp-med)*(par-med))/((perp + med)*(par*par -2.*par*med + med*med))
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
    #multiplies aiz by 1 = (par-med)/(par-med) to get zero instead of
    #divergence

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

T,Y = np.meshgrid(t,y)

U,Y0 = np.meshgrid(u,y0)
		
# Calculate N=0 term        
f_term0 = (1./(Y0*np.sqrt(1.0 - (1./Y0)))) *U*U*U*U\
        *np.exp(-2.* Y0 *U)\
        *(2.*(1. + 3.*a[0])*(1. + 3.*a[0])\
        +(1.-a[0])*(1.-a[0]))
Ft_term0  = np.trapz(f_term0, U, axis = 1)
Fty_term0 = np.trapz(Ft_term0,y0)
A_term0   = (1./2) *kbT*(1./(12.*np.pi))* delta[0]*delta[0]*Fty_term0

#print 'Ft_term0 = ',Ft_term0
#print 'Fty_term0 = ',Fty_term0
# Calculate 1 < N < 500 terms

A = np.zeros(shape = (len(Ls),len(ns)))
for i,L in enumerate(Ls):
    sum_A = np.zeros(shape = (len(Ls),len(ns)))
    for j,n in enumerate(ns):
        #print 'N,L = ',(n,L)
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A
        f = (1./(Y*np.sqrt(1.0 - (1./Y))))\
                * np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1.))\
                * (T / np.sqrt(T*T + 1.))\
                * (2. * (1. + 3.*a[j]) * (1. + 3.*a[j]) * T*T*T*T\
                + 4. * (1. + 2.*a[j] + 2.*a[j] + 3.*a[j] * a[j]) *T*T\
                + 4.*(1.+a[j])*(1.+a[j])\
                + (T*T*T*T + 4.*T*T + 4.)*(1.-a[j])*(1.-a[j]))
        #print 'f = ',f
        #f = (1./np.sqrt(Y*Y -1.0))\
        #        * np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1.))\
        #        * (T / np.sqrt(T*T + 1.))\
        #        * (2. * (1. + 3.*a[j]) * (1. + 3.*a[j]) * T*T*T*T\
        #        + 4. * (1. + 2.*a[j] + 2.*a[j] + 3.*a[j] * a[j]) *T*T\
        #        + 4.*(1.+a[j])*(1.+a[j])\
        #        + (T*T*T*T + 4.*T*T + 4.)*(1.-a[j])*(1.-a[j]))
       
        #Ft = np.trapz(f , T, axis = 1)
        Ft = simps(f , T, axis = 1)
        #print 'Ft = ', Ft
        #Ft = np.sum(f , axis = 1)
        #Fty = np.sum(Ft)
        #Fty = np.trapz(Ft,y)
        Fty = simps(Ft,y)
        #print 'Fty = ',Fty
        #print 'inner L = ',i
        #print 'inner j,N = ',j,n
        # Matsubara terms 
        A[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Fty
    #print 'outer L = ',i
    #print 'outer j,N = ',j,n
    #sum_A = (kbT/(12.* np.pi)) * np.sum(A, axis = 1)
    #print '*********************************', np.shape(sum_A),sum_A
sum_A = (kbT/(12.* np.pi)) * np.sum(A, axis = 1)
#print 'f_term0 = ',f_term0
#print 'Fty_term0 = ',Fty_term0
#print 'A_term0 = ', (kbT/(12.*np.pi))*A_term0
#print 'Ft = ',Ft
#print 'Fty = ',Fty
print 'sum_A = ',sum_A

np.savetxt('65w65_prw_sum_A_mvsum.txt',sum_A)
np.savetxt('65w65_prw_Ls_mvsum.txt',Ls)

# Plots
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
pl.loglog(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_par)
pl.title(title('6','5','parallel'))
pl.legend(loc = 'best')
#pl.savefig(svfig('9','1','loglog_parallel'))
pl.show()
#pl.figure()
#pl.loglog(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_par)
#pl.title(title('6','5','parallel'))
#pl.legend(loc = 'best')
#pl.savefig(svfig('6','5','loglog_parallel'))
#pl.show()
#
#pl.figure()
#pl.semilogy(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.title(title('6','5','parallel'))
#pl.legend(loc = 'best')
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.savefig(svfig('6','5','parallel'))
#pl.show()
#pl.legend(loc = 'best')
#pl.show()
#
