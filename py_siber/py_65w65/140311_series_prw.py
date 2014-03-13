#!/usr/bin/python
import numpy as np
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

eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output_eV.txt')
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # parallel,axial
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # water as intervening medium
#eiz_w[0] = eiz_w[1] #NOTE: there is a jump from first val down to second val

x00_A0, y00_A0 = np.loadtxt('data/65-W-65-PARA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x00_A2, y00_A2 = np.loadtxt('data/65-W-65-PARA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium

r_1 = 1.0e-9
r_2 = 1.0e-9
c = 2.99e8 # in m/s
Temp = 297. 
kb  = 1.3807e-23 # in J/K
coeff = 2.411e14 # in rad/s
# NOTES:
# at RT, 1 kT = 4.11e-21 J
# 1 eV = 1.602e-19 J = 0.016 zJ
# h_bar_eV = 6.5821e-16 eVs
# h_bar = 1. #1.0546e-34 #in Js
#kb = 8.6173e-5  # in eV/K
# z_n_eV = (2*pi*kT/h_bar)n 
#	= (0.159 eV) / (6.5821e-16 eVs) 
#	= n*2.411e14 rad/s
# z_n_J = (2*pi*kT/h_bar)n 
#	= (1.3807e-23 J/K) / (1.0546e-34 Js))*n
#	= n*2.411e14 rad/s
#coeff = 0.159 # in eV w/o 1/h_bar

ns = np.arange(0.,500.)
z = ns * coeff
ls = np.linspace(1.0e-9, 1.0e-6, 200)
#thetas = np.linspace((0.0001)*np.pi,(1./2)*np.pi,25)
dt = 1.0
ts = np.arange(1.0,201.,dt)
dh = 1.0
hs = np.arange(1.01,201.,dh)

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def g0(a,eizw,L, N,yy,time):
	term0 = ( time    / np.sqrt(time*time+1.)   )
	term1 =	( time**(4.) * 2.0*(1. + 3.*a)*(1.+3.*a)     )
	term2 = ( time**(2.) * 4.0*(1. + 2.*a+2.*a+3.0*a*a))
	term3 = (          4.0*(1. + a)*(1. + a)         )
	#term4 = (-2.0 * yy * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	term4 = yy * np.sqrt(eizw) * L * L * N * (coeff/c) * np.sqrt(time*time + 1.0)	#NOTE where the extra L comes from subst of y to y/l
	#term5 = (1./(np.sqrt(yy*yy - 1)))
	term5 = (1./(np.sqrt(yy*yy - 1.)))
#	print 'ys term0', term0
#	print 'ys term1', term1
#	print 'ys term2', term2
#	print 'ys term3', term3
	print 'ys term4', np.exp(term4)
	print 'ys term4', term4
#	print 'ys term5', term5
	#print '----'
	#return (term0) * np.exp(term4)*( (term1) + (term2) + (term3)) * term5
	return term0 * (term1 + term2 + term3) * term5 * ( 1.0 - 2.0*term4 + 2.0* term4*term4 - (4./3)*term4*term4*term4 + (2./3)*term4*term4*term4*term4)  

def g2(a,eizw, L, N, yy,time): 
	term0 = (time    / np.sqrt(time*time+1.0)   ) 
	term1 = ((1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0))                     
	term2 = yy * np.sqrt(eizw) * L * L * N * (coeff/c) * np.sqrt(time*time + 1.0)	
	#term2 = (-2.0 * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	term3 = (1./(np.sqrt(yy*yy - 1.)))
	#term2 = (-2.0 * np.sqrt(eizw) * L * coeff * N * np.sqrt(time*time + 1.0) * (1./c)) 
	#term2 = (np.exp(-2.0 * yy * np.sqrt(eizw)* L * N * (1./c) * np.sqrt(time*time + 1.0)))** (coeff) 
	#term2 = yy * np.sqrt(eizw) * L * N * (coeff/c) * np.sqrt(time*time + 1.0)	
	#term3 = (1./(np.sqrt(yy*yy - 1)))
#	print 'y_2s term0', term0
#	print 'y_2s term1', term1
#	#print 'y_2s term2', np.exp(term2)
#	print 'y_2s term2', term2
#	print 'y_2s term3', term3
#	print '----'
	#return term0 * term1* np.exp(term2) * term3
	return term0 * term1 * term3 * ( 1.0 - 2.0*term2 + 2.0* term2*term2 - (4./3)*term2*term2*term2 + (2./3)*term2*term2*term2*term2)  

def As(eizz,eizw,L,N,Y): 
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * ((np.sqrt(eizw))**(5.)) * (L**(5.)) * (coeff**(5.)) * (N**(5.)) / (c**(5.))) 
	#term2 = np.exp(Y * eizw **(5./2) *((L**5) / (c**5)) * (coeff*N)**5 ) 
	#term3 = (L**5) / (c**5)
	#print 'As term1 = ', term1
#	print 'As term2 = ', term2
	#print '----'
	#return  np.log(term1 +  term2)# * term3
	return  term1 *  term2# * term3

def A_2s(eizz,eizw, L , N ,Y):		
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * ((np.sqrt(eizw))**(5.)) * (L**(5.)) * (coeff**(5.)) * (N**(5.)) / (c**(5.))) 
	#term3 = (L**5) / (c**5)
	#print 'A_2s term1 = ', term1
#	print 'A_2s term2 = ', term2
	#print '----'
	#return np.exp(term1 + term2)# * term3
	return term1 * term2# * term3

aiz = np.zeros(len(ns))
aiz = Aiz(eiz_x,eiz_z, eiz_w) 

H,T  = np.meshgrid(hs,ts)
#g0 = np.zeros(shape=(len(hs),len(ts)))
#g2 = np.zeros(shape=(len(hs),len(ts)))
G0dt = np.zeros(len(ts))
G2dt = np.zeros(len(ts))

G0 = np.zeros(shape=(len(ns),len(ls)))
G2 = np.zeros(shape=(len(ns),len(ls)))

A   = np.zeros(shape=(len(ns),len(ls)))
A_2 = np.zeros(shape=(len(ns),len(ls)))

EL = np.zeros(len(ls))
G_l_t_dt = np.zeros(len(ls))

for k,length in enumerate(ls):
    	#print "on l[k], k =%d of %d"%(k,len(ls))
	sum_A = np.empty(len(ls))
	sum_A_2 = np.empty(len(ls))
	for j,n in enumerate(ns):
		#print 'dt dh Integrand g0 = ',j, g0(aiz[j],eiz_w[j],length,n, H, T)
		#print 'dt dh Integrand g2 = ',j, g2(aiz[j],eiz_w[j],length,n, H, T)

		# Integral:
		g0dt = g0(aiz[j],eiz_w[j],length, n, H, T)
		G0dt = trapz(g0dt, ts, axis = 0)
		G0   = trapz(G0dt, hs)

		g2dt = g2(aiz[j],eiz_w[j],length,n, H, T)
		G2dt = trapz(g2dt, ts, axis = 0)#(aiz[j],eiz_w[j],length,n, H, T),T)
		G2   = trapz(G2dt, hs, axis = 0)

		#print 'dt Integral G0 = ',j, G0
		#print 'dt Integral G2 = ',j, G2
		#print '----'
		#print 'N terms for A0 = '  , As(eiz_z[j],eiz_w[j],length,n,G0)
		#print 'N terms for A2 = ', A_2s(eiz_z[j],eiz_w[j],length,n,G2)
		#print '----'

		A[j,k]   =   As(eiz_z[j],eiz_w[j],length,n,G0)
		A_2[j,k] = A_2s(eiz_z[j],eiz_w[j],length,n,G2)# * np.cos(2.0*theta)  		

	A[0,k] = (1./2)*A[0,k]      
	A_2[0,k] = (1./2)*A_2[0,k]  

	sum_A = np.sum(A,axis=0)

	#print 'sum of A0 = ', k,j,sum_A
	sum_A_2 = np.sum(A_2,axis=0)
	#print 'sum of A2 = ', k,j,sum_A_2
	#print '----'
	#print 'shape sum_A_2 = ', np.shape(sum_A_2)
	#sys.exit()
for k,length in enumerate(ls):
	EL[k] = 1./(length*length*length*length*length)
	G_l_t_dt[k] = (1.602e-19 / 4.11e-21) * (1./32) * EL[k]*r_1*r_1*r_2*r_2*(1./(12*np.pi))*(sum_A[k] + sum_A_2[k])

np.savetxt('g_prw.txt',G_l_t_dt)
#pl.figure()
#pl.plot(ns,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
#pl.plot(ns,eiz_y, color = 'g', label = r'$\varepsilon_{\hat{y}}(i\zeta_{N})$')
#pl.plot(ns,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{N})$')
##pl.plot(ns,eiz_w, color = 'c', label = r'$\varepsilon_{vac}(i\zeta_{N})$')
#pl.plot(ns,eiz_w, color = 'c', label = r'$\varepsilon_{water}(i\zeta_{N})$')
#pl.xlabel(r'$N$', size = 20)
#pl.ylabel(r'$\varepsilon(i\zeta)$', size = 20)
#pl.legend(loc = 'best')
#pl.title(r'$\mathrm{CG-10\, DNA}$', size = 20)
#pl.axis([0,500,0.9,2.6])
#pl.savefig('plots/par_ret_water/eiz.pdf' )
#show()

pl.figure()
# convert eV to zJ with 1.602e2 zJ/eV and here J = 1e9*N*m = N*nm
#pl.loglog(ls,(1.0e9*1.602e2*kb*Temp/(12.*np.pi))*(sum_A + sum_A_2),'b', label = r'$\mathcal{A^{(0)}+A^{(2)}}$')
#pl.loglog(ls,(1.0e9*1.602e2*kb*Temp/(12.*np.pi))*sum_A            ,'g', label = r'$\mathcal{A^{(0)}}(\ell,\theta = 0)$')
#pl.loglog(ls,(1.0e9*1.602e2*kb*Temp/(12.*np.pi))*sum_A_2          ,'r', label = r'$\mathcal{A^{(2)}}(\ell,\theta = 0)$')


pl.loglog(1.e9*ls,(1.e21*kb*Temp/(12.*np.pi))*(sum_A + sum_A_2),'b', label = r'$\mathcal{A^{(0)}+A^{(2)}}$')
pl.loglog(1.e9*ls,(1.e21*kb*Temp/(12.*np.pi))*sum_A            ,'g', label = r'$\mathcal{A^{(0)}}(\ell,\theta = 0)$')
pl.loglog(1.e9*ls,(1.e21*kb*Temp/(12.*np.pi))*sum_A_2          ,'r', label = r'$\mathcal{A^{(2)}}(\ell,\theta = 0)$')

pl.loglog(x00_A0, y00_A0,'k-' , label = r'GH $\mathcal{A^{(0)}}(\ell,\theta = 0)$')
pl.loglog(x00_A2, y00_A2,'k--', label = r'GH $\mathcal{A^{(2)}}(\ell,\theta = 0)$')

pl.xlabel(r'$\ell\,\,\,\rm{[m]}$', size = 24)
pl.ylabel(r'$\mathcal{A(\ell)}$', size = 20)
#pl.title(r'$\mathrm{Hamaker \, coeff.s \,:\,skewed,\,retarded,\,water}$', size = 20)
pl.legend(loc = 'upper right')
#pl.axis([1e-9,1e-6,1e-24,1e-19])
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/par_ret_water/140306_65w65_par_ret_A.pdf')
show()

pl.figure()
pl.loglog(ls, G_l_t_dt)#, label = labels[i])
pl.xlabel(r'$\ell\,\,\mathrm{[m]}$', size = 24)
pl.ylabel(r'$-g(\ell,\theta = 0)\,\,\mathrm{[k_{B}T/nm]}$', size = 20)
#pl.axis([1.9e-9, 8.1e-9,1.5e7,1.5e10])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,parallel,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/par_ret_water/140306_65w65_par_ret_g_vs_l.pdf')
show()



