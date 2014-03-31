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
eiz_w[0] = eiz_w[1] #NOTE: there is a jump from first val down to second val

r_1 = 1.0e-9
r_2 = 1.0e-9
c = 2.99e8 # in m/s
Temp = 297 
kb  = 1.3807e-23 # in J/K
coeff = 2.411e14 # in rad/s

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
	term0 = ( time    / np.sqrt(time*time+1.0)               )
	term1 =	( time**4 * 2.0*(1. + 3.*a)*(1.+3.*a)     )
	term2 = ( time**2 * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a))
	term3 = (          4.0*(1. + a)*(1.0 + a)         )
	term4 = (-2.0 * yy * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	term5 = (yy*yy - 1)**(-1./2)
	return (term0) * np.exp(term4)*( (term1) + (term2) + (term3)) * term5

def g2(a,eizw, L, N, yy,time): 
	term0 = (time    / np.sqrt(time*time+1.0)                               ) 
	term1 = ((1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0))                     
	term2 = (-2.0 * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	term3 = (yy*yy - 1)**(-1./2)
	return term0 * term1* np.exp(term2) * term3

def As(eizz,eizw,L,N,Y): 
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * eizw **(5./2) * (coeff*N)**5 * L**5 / (c**5))
	return  term1 *  term2

def A_2s(eizz,eizw, L , N ,Y):		
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * eizw *eizw*np.sqrt(eizw) * (coeff*N)**5 * L**5 / (c**5))
	return (term1 * term2)

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

	sum_A = np.empty(len(ls))
	sum_A_2 = np.empty(len(ls))

	for j,n in enumerate(ns):
	

#	for h in hs:
	#		term1 = 1./np.sqrt(h**2-1)
	#		for t in ts:
	#			intermediate_integral(t,h) += np.exp(-2.0 * h * np.sqrt(eiz_w)* length * coeff * n / c * np.sqrt(t*t + 1.0))*((t/ np.sqrt(t*t+1.0)) + ( t**4 * 2.0*(1. + 3.*a)*(1.+3.*a)) + ( t**2 * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a)) + (4.0*(1. + a)*(1.0 + a)))

	#		h_integrand += term1*intermediate_integral()

	#			g0dt = g0(aiz[j],eiz_w[j],length,n,h, t)
	#
	#			g2dt = g2(aiz[j],eiz_w[j],length,n, H, T)

	#	## Integral:
	#	#g0dt = g0(aiz[j],eiz_w[j],length,n, H, T)
	#	#G0dt = trapz(g0dt, x=ts, axis = 1)
	#	#G0   = trapz(G0dt, x=hs)

	#	#g2dt = g2(aiz[j],eiz_w[j],length,n, H, T)
	#	#G2dt = trapz(g2dt, x=ts, axis = 1)#(aiz[j],eiz_w[j],length,n, H, T),T)
	#	#G2   = trapz(G2dt, x=hs)
	#	#print 'dt Integral   G0 = ',i,k,j, G0
	#	#print 'dt Integral G2 = ',i,k,j, G2
	#	#print '----'
	#	#print 'N terms for A0 = '  , As(eiz_z[j],eiz_w[j],length,n,G0)
	#	#print 'N terms for A2 = ', A_2s(eiz_z[j],eiz_w[j],length,n,y_2)
	#	#print '----'

		A[j,k]   =   As(eiz_z[j],eiz_w[j],length,n,G0)
		A_2[j,k] = A_2s(eiz_z[j],eiz_w[j],length,n,G2)# * np.cos(2.0*theta)  		

	A[0,k]   *= 0.5  
	A_2[0,k] *= 0.5

	sum_A = np.sum(A,axis=0)

	#print 'sum of A0 = ', k,j,sum_A
	sum_A_2 = np.sum(A_2,axis=0)
	#print 'sum of A2 = ', k,j,sum_A_2
	#print '----'
	#print 'shape sum_A_2 = ', np.shape(sum_A_2)
	#sys.exit()
for k,length in enumerate(ls):
	EL[k] = 1./(length*length*length*length*length)
	#G_l_t_dt[k] = (1.602e-19 / 4.11e-21) * (3.*np.pi/8) * EL[k]*r_1*r_1*r_2*r_2*(1./(12*np.pi))*(sum_A[k] + sum_A_2[k])
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
pl.loglog(ls,(kb*Temp/(12.*np.pi))*(sum_A + sum_A_2))#,'b-', label = r'$\mathcal{A^{(0)}}(\ell)$')
#pl.loglog(ls,(kb*Temp/(12.*np.pi))*sum_A_2,'g-', label = r'$\mathcal{A^{(2)}}(\ell)$')
#pl.loglog(ls,(kb*T/32)*A2_theta,':', label = r'$\mathcal{A^{(2)}}(\ell)cos(2\theta)$')
pl.xlabel(r'$\ell\,\,\,\rm{[m]}$', size = 24)
pl.ylabel(r'$\mathcal{A(\ell)}$', size = 20)
#pl.title(r'$\mathrm{Hamaker \, coeff.s \,:\,skewed,\,retarded,\,water}$', size = 20)
pl.legend(loc = 'upper right')
#pl.axis([1e-9,1e-6,1e-24,1e-19])
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/par_ret_water/140306_93w93_par_ret_A.pdf')
pl.show()

#pl.figure()
#pl.loglog(ls,(kb*T/32)*sum_A,'b-', label = r'$\mathcal{A^{(0)}}$')
#pl.loglog(ls,(kb*T/32)*sum_A_2,'b--', label = r'$\mathcal{A^{(2)}}$')
#pl.xlabel(r'$\mathrm{separation}\,\ell\,\,\,\rm{[m]}$', size = 20)
#pl.ylabel(r'$\mathrm{\mathcal{-A^{(0)},\,\,-A^{(2)}}}$', size = 20)
##pl.title(r'$\mathrm{Hamaker \, coeff.s \,:\,parallel,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'lower right')
##pl.axis([2e-9,1e-8,1e-24,1e-18])
#pl.savefig('plots/par_ret_water/140306_290w290_par_A0_A2.pdf')
#show()

#pl.figure()
##pl.loglog(thetas, G_l_t_dt)#, label = labels_l[k])
#pl.semilogy(thetas, G_l_t_dt)
#pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
#pl.ylabel(r'$G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([(1./25)*np.pi,(3./4)*np.pi,105,135])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,parallel,\,retarded,\,water}$', size = 20)
##pl.legend(loc = 'best')
#pl.savefig('plots/par_ret_water/semilog_G_vs_theta.pdf')
#show()

pl.figure()
pl.loglog(ls, G_l_t_dt)#, label = labels[i])
pl.xlabel(r'$\ell\,\,\mathrm{[m]}$', size = 24)
pl.ylabel(r'$-g(\ell,\theta = 0)\,\,\mathrm{[k_{B}T]}$', size = 20)
#pl.axis([1.9e-9, 8.1e-9,1.5e7,1.5e10])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,parallel,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/par_ret_water/140306_93w93_par_ret_g_vs_l.pdf')
pl.show()


