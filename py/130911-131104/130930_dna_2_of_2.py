#!/usr/bin/python

import numpy as np
from numpy import log,exp
#from scipy import *
import scipy.optimize as opt
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import axis as ax
import matplotlib               
from pylab import *
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

x_t ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

eiz_x = np.loadtxt('data/eiz_x_output.txt')
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt')

# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV #(2.41*1e14) # in rad/s
ns = arange(0,250)
z = ns * coeff
kb = 8*1e-5 #in eV/K
T = 300 # r.t. in K
c = 3.0e8
r_1 = 1e-9
r_2 = 1e-9
t = np.linspace(0,100,100)
G_const = -1.0*r_1*r_2
#ls = 8e-10#[1e-6 , 1e-5 , 1e-4]
ls = linspace(8e-10 , 1e-9 , 10)
thetas = linspace(np.pi/10, 2*np.pi, 3)#/4, np.pi/3, np.pi]
eiz_m = 1.0

#height= np.zeros(shape=(len(ls),len(thetas)))
aiz = np.zeros(len(eiz_y))
def a(perp, par):
    return 2*(perp-1)/((perp+1)*(par-1))
aiz = a(eiz_y,eiz_x)

pl.figure()
pl.plot(ns,aiz)
pl.show()





#def a(perp, par):
#	""" Implicit definition of a(perp, par, med) """
#	return 2*(perp - 1)*1 / ((perp + 1)(par - 1))
#aiz = [a(eizy,eizx) for eizy,eizx in zip(eiz_y,eiz_x)]# 2.0*(eiz_y[m] - 1.0)*1.0 / ((eiz_y[m] + 1.0)(eiz_x[p] - 1.0))# for eizy,eizx in zip(eiz_y,eiz_x)]

def g_0(thetas,aa,time):
	""" Implicit definition of a(perp, par, med) """
	return 2.0*(1.0/sin(thetas)) * ((1. +3.*aa)*(1. + 3.*aa)*time*time*time*time + 2.*(1. + 2.0*aa + 2.0*aa + 3.0*aa*aa) *time*time + 2.0 * (1.0+aa) * (1.0 + aa))

def g_2(thetas,aa,time):
	""" Implicit definition of a(perp, par, med) """
	return 2.0*(cos(2.0*thetas)/sin(thetas)) * ((1.0 - aa)*(1.0 - aa)*(time*time + 2.0)*(time*time + 2.0))

#def f(p, phi, phib, df):
#	""" Implicit definition of P(phi, dF) """
#	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
#		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))
#
#def P(phi, phib, df):
#	""" Numerically solve for partition coefficient as a
#	    function of \phi_s """
#	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
#		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
#	else:
#		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson
#
lts = [':','--','-']
markers = ['x', '^', 'd']
for j in range(len(ns)):
    sum_Li2[j] = 0.0
    for k in range(len(ls)):
    	for k in range(len(ls)):
		g0 = g_0(theta,aiz,t) 
		g2 = g_2(theta,aiz,t)
		int_g_arg = [t/(t*t+1)* np.exp(-2*1*(t*t+1)**(1./2)*coeff*n*l/c)*(g_0(theta,a_iz,t) + g_2(theta,a_iz,t)) for n,a_iz in zip(ns,aiz)]
		int_g = trapz(int_g_arg,t)	
		pl.plot(ns,int_g)
		A = (kb)*T/32*((eiz_x - 1)/1)*((eiz_x - 1)/1)*1*((coeff*n)*(coeff*n)*l*l/c/c)*((coeff*n)*(coeff*n)*l*l/c/c)*int_g


	        sum_Li2[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**m)/(m**2)
	    dF[j] = -sum_Li2[j]\
	    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
	    *(eiz_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
	pl.figure()
for i,n in enumerate(ns):
	for j,l in enumerate(ls): 
		for k,theta in enumerate(thetas):
	#	for q,t in enumerate(ts):
			g0 = g_0(theta,aiz,t) 
			g2 = g_2(theta,aiz,t)
			int_g_arg = [t/(t*t+1)* np.exp(-2*1*(t*t+1)**(1./2)*coeff*n*l/c)*(g_0(theta,a_iz,t) + g_2(theta,a_iz,t)) for n,a_iz in zip(ns,aiz)]
			int_g = trapz(int_g_arg,t)	
		pl.plot(ns,int_g)
		A = (kb)*T/32*((eiz_x - 1)/1)*((eiz_x - 1)/1)*1*((coeff*n)*(coeff*n)*l*l/c/c)*((coeff*n)*(coeff*n)*l*l/c/c)*int_g
#	A_sum = sum(A)
		pl.plot(ls,A)
       # 	height[j,k] = A
#    	        except: continue
#    	    	if j==0:
#    	    		labels=r'$\phi_{b} = %.2f$' % phi_b
#    	    	else:
#    	    		labels=None
#    	    	if i==2:
#    	    		labelss=r'$\Delta f = %.1f$' % df
#    	    	else:
#    	    		labelss=None
#    	    	print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
#    	    	op = osmotic_pressure(phi,phi_b)

######	X,Y = meshgrid(ls,thetas)

#####figure()
#contourf(X,Y,height,100); clim(0,4)#,cmap = hot())
#levels = [0.1,0.5,1.0,1.2,1.3,1.4,1.5,2.0]#,0.80,1.00]
######contour(X,Y,height)
#CS = contour(X,Y,height,levels, colors = 'k')
#man_loc = [(1,1),(2,2),(3,3),(4,4)]
#clabel(CS, inline =1,fmt = '%2.1f', fontsize = 18,color = 'k', manual = man_loc)
#contourf(X,Y,height, 1000); clim(0,4)#;cbar = colorbar();
#CS=contour(X,Y,height, levels)#, colors = 'k')
#clabel(CS,fmt = '%2.1f',colors = 'k',fontsize = 18)
#xlabel(r'$\mathrm{\Phi_{PEG\,400}\,\,added\, to\, fixed\,\,\Phi_{PEG\,3500}=\,0.20}$', size = 18)#\,=\,\frac{2}{\pi}C$', size = 20)
#ylabel(r'$\Delta f$', size = '20')#\,\,\,Energy cost to enter pore$', size = 20)
#cbar.ax.set_ylabel(r'$P\,=\,\frac{\Phi_{400}^{in}}{\Phi_{400}^{out}}$', size = 20)
#cbar.add_lines(CS)
#axes([0,0.78,0,3])

#####savefig('plots/130930_test')
#cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
#cbar.add_lines(CS)
#####show()
####close()
