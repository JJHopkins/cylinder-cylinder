#!/usr/bin/env python
import numpy as np
from scipy.integrate import romb
from scipy.integrate import simps #quad
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

#Ls = np.arange(0.,1e2)#1e-9,101e-9,10e-9)  # separation distance between 2 cyclinders
Ls = np.arange(1e-9,1e-6,1e-9)  # separation distance between 2 cyclinders
#t  = np.logspace(0.,200,201)#,100)
#u  = np.logspace(1.,200,201)#,99) 
#y  = np.logspace((1.+1e-5),200,200)#,98)
#y0 = np.logspace((1.+1e-5),200,200)#,98)

t  = np.linspace(0.,2.**10, 1.+2.**10)#200,201)#,100)
#u  = np.linspace(0.,200,201)#,99) 
u  = np.linspace(0.,2.**10, 1.+2.**10)#,99) 
y  = np.linspace((1.+1e-5),2.**10, 1.+2.**10)#200,200)#,98)
y0 = np.linspace((1.+1e-5),2.**10,1.+2.**10)#,98)
#y0 = np.linspace((1.+1e-5),200,200)#,98)

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
A = np.zeros(shape = (len(Ls),len(ns)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

T,Y = np.meshgrid(t,y)

U,Y0 = np.meshgrid(u,y0)
		
# Calculate N=0 term        
f_term0 = (1./(Y0*np.sqrt(1.0 - (1./Y0)))) *U*U*U*U\
        *np.exp(-2.* Y0 *U)\
        *(2.*(1. + 3.*a[0])*(1. + 3.*a[0])\
        +(1.-a[0])*(1.-a[0]))
#Ft_term0 = simps(f_term0, U , even = 'last')
Ft_term0 = romb(f_term0, axis = 1)#U , even = 'last')
#Ft_term0 = np.sum(f_term0, axis = 1)
#Ft_term0 = np.trapz(f_term0, dx = 2,axis = 1)
#Ft_term0 = np.trapz(f_term0, u)# axis = 1)
#Fty_term0 = np.sum(Ft_term0)
#Fty_term0 = np.trapz(Ft_term0)
#Fty_term0 = simps(Ft_term0, y0, even = 'last')
Fty_term0 = romb(Ft_term0, axis = 0)#, y0, even = 'last')
#Fty_term0 = np.trapz(Ft_term0,y)
#print 'Ft_term0 = ',Ft_term0
#print 'Fty_term0 = ',Fty_term0
# Calculate 1 < N < 500 terms
for i,L in enumerate(Ls):
    #print L
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A
        f = (1./(Y*np.sqrt(1.0 - (1./Y))))\
                * np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1.))\
                * (T / np.sqrt(T*T + 1.))\
                * (2. * (1. + 3.*a[j]) * (1. + 3.*a[j]) * T*T*T*T\
                + 4. * (1. + 2.*a[j] + 2.*a[j] + 3.*a[j] * a[j]) *T*T\
                + 4.*(1.+a[j])*(1.+a[j])\
                + (T*T*T*T + 4.*T*T + 4.)*(1.-a[j])*(1.-a[j]))
       
        #f = (1./np.sqrt(Y*Y -1.0))\
        #        * np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1.))\
        #        * (T / np.sqrt(T*T + 1.))\
        #        * (2. * (1. + 3.*a[j]) * (1. + 3.*a[j]) * T*T*T*T\
        #        + 4. * (1. + 2.*a[j] + 2.*a[j] + 3.*a[j] * a[j]) *T*T\
        #        + 4.*(1.+a[j])*(1.+a[j])\
        #        + (T*T*T*T + 4.*T*T + 4.)*(1.-a[j])*(1.-a[j]))
       
        #Ft = simps(f , T)
        #Ft = np.trapz(f , T)
        #Ft = np.sum(f , axis = 1)
        Ft = romb(f , axis = 1)#U , even = 'last')
        Fty =romb(Ft, axis = 0)#U , even = 'last')
        #Fty = np.sum(Ft)
        #Fty = simps(Ft,y)
        #Fty = np.trapz(Ft,y)
        # Matsubara terms 
    A[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Fty
    A[i,0] = (1./2) * delta[0]*delta[0]*Fty_term0
sum_A = (kbT/(12.* np.pi)) * np.sum(A, axis = 1)
#print 'f_term0 = ',f_term0
#print 'Fty_term0 = ',Fty_term0
print 'A_term0 = ', (kbT/(12.*np.pi))*A[:,0]
#print 'Ft = ',Ft
#print 'Fty = ',Fty
print 'sum_A = ',sum_A

np.savetxt('65w65_prw_sum_A.txt',sum_A)
np.savetxt('65w65_prw_Ls.txt',Ls)

#pl.figure()
#pl.loglog(Ls,sum_A,label=r'$\mathcal{A}(\ell= %1.1f nm)= %3.2f \,\,\, \mathcal{A}(\ell= %1.1f nm)= %3.2f $'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[3],1e21*sum_A[3]))
#pl.legend(loc = 'best')
#pl.show()

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

#pl.figure()
#pl.plot(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_par)
#pl.title(title('9','1','parallel'))
#pl.legend(loc = 'best')
#pl.savefig(svfig('9','1','loglog_parallel'))
pl.show()
pl.figure()
pl.loglog(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_par)
pl.title(title('6','5','parallel'))
pl.legend(loc = 'best')
pl.savefig(svfig('6','5','loglog_parallel'))
pl.show()
#
pl.figure()
pl.semilogy(1e9*Ls,1e21*sum_A,label=r'$\mathcal{A}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A[0],1e9*Ls[1],1e21*sum_A[1]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.title(title('6','5','parallel'))
pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig(svfig('6','5','parallel'))
pl.show()
pl.legend(loc = 'best')
pl.show()
#

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




