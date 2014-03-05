#!/usr/bin/python
import numpy as np
from numpy import log,exp
#from scipy import *
import scipy.optimize as opt
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import axis as ax
import matplotlib               
from pylab import *
#from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from mpl_toolkits.mplot3d import Axes3D

x_t_eV ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x_eV ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y_eV ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z_eV ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

x_t =x_t_eV /6.582e-16
x_x =x_x_eV /6.582e-16
x_y =x_y_eV /6.582e-16
x_z =x_z_eV /6.582e-16

# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 2.41e14 # in rad/s
ns = arange(0.0,500.0)
T = 300.0 # r.t. in K
kb_J =1.3806488e-2# in zJ/K # 1.3806488e-23 # in J/K
z = ns * coeff
r_1 = 0.5e-9
r_2 = 0.5e-9
#r_1 = 0.5
#r_2 = 0.5
t = np.linspace(0.0,10000.0,100)
ls = [4.0e-9]#linspace(0.5e-9, 20e-9, 100)
c = 2.99e8
#thetas = [np.pi/7.0, np.pi/5.0, np.pi/4.0, np.pi/3.0, np.pi/2.00001]#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
thetas = linspace(np.pi/100.0,np.pi/2.0,50)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
#thetas = linspace(np.pi/3.0,np.pi/2.0,4)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
#thetas = [0.000001, np.pi/2.000001]#,np.pi/4.0,np.pi/2.0]#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#

eiz_x = np.loadtxt('data/eiz_x_output.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') # parallel,axial
eiz_m = 1.0

aiz = np.zeros(len(eiz_z))
def a(perp, par):
    return 2.0*(perp-1.0)/((perp+1.0)*(par-1.0))
aiz = a(eiz_x,eiz_z)

lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']

int_g_arg = np.empty(len(ns))
int_g = np.empty(len(ns))
int_g_2 = np.empty(len(ns))
A = np.empty(len(ns))
A_2 = np.empty(len(ns))
A0 = np.empty(len(ns))
A2 = np.empty(len(ns))
sum_A = np.empty(len(ls))
Sum_A = np.empty(len(ls))
sum_A_thetas = np.empty(len(thetas))
sum_A_2 = np.empty(len(ls))
Sum_A_2 = np.empty(len(ls))
sum_A_2_thetas = np.empty(len(thetas))
G = np.empty(len(ls))
G_2 = np.empty(len(ls))
added = np.empty(len(ls))
added_thetas = np.empty(len(thetas))

pl.figure()
for p in range(len(thetas)):
	for k in range(len(ls)):
		for j in range(len(ns)):
			int_g[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(2.0*((1.+3.*aiz[j])*(1.+3.*aiz[j])*t*t*t*t+2.*(1.+2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])*t*t+2.0*(1.0+aiz[j])*(1.0 + aiz[j]))),t)

		#	y = np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)

			int_g_2[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
			(((1.0 - aiz[j])*(1.0 - aiz[j])*(t*t + 2.0)*(t*t + 2.0))),t)

		#	A[j] = (kb)*T/32.0*((eiz_x[j]-1.0)/1.0)*((eiz_x[j]-1.0)/1.0)*\
		#	1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
		#	int_g[j]

		#	A_2[j] = (kb)*T/32.0*((eiz_x[j]-1.0)/1.0)*\
		#	((eiz_x[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
		#	int_g_2[j]	
			
			A[j] = kb_J*T/32.0*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*\
			1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g[j]

			A_2[j] = kb_J*T/32.0*((eiz_z[j]-1.0)/1.0)*\
			((eiz_z[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g_2[j]	



			Sum_A[k] += kb_J*T/32.0*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*\
			1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g[j]

			Sum_A_2[k] += kb_J*T/32.0*((eiz_z[j]-1.0)/1.0)*\
			((eiz_z[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g_2[j]	
			
			#sum_A[k] = sum(A[j])
			sum_A_thetas[p] +=  kb_J*T/32.0*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*\
			1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g[j]
			
			#sum_A_2[k] = sum(A_2[j])
			sum_A_2_thetas[p] += kb_J*T/32.0*((eiz_z[j]-1.0)/1.0)*\
			((eiz_z[j]-1.0)/1.0)*1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/c/c)*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
			int_g_2[j]	

		added_thetas[p] = -(-((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_thetas[p])+\
		(-((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/\
		(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*sum_A_2_thetas[p])	

		G = -((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A

		G_2 = -((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A_2

		added = (-((np.pi*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A)+(-((np.pi*np.pi*r_1*r_1*r_2*r_2*cos(2.0*thetas[p]))/\
		(2.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A_2)	

pl.plot(thetas,added_thetas)#, color = colors[p])#, linestyle = lts[k])	

#pl.text(0.1,-0.15e-8,r'$\theta = %.2f \, radians$' % thetas[0], color='r')
#pl.text(0.1,-0.30e-8,r'$\theta = %.2f \, radians$' % thetas[1], color='g')
#	pl.text(0.1,-0.45e-8,r'$\theta = %.2f \, radians$' % thetas[2], color='b')
pl.legend()
pl.xlabel(r'$\theta\,\,\,[radians]$', size = 24)
pl.ylabel(r'$G(l,\theta)$', size = 24)
#pl.legend()
pl.title(r'Eqn 12 $G(l,\theta)\,with\,\varepsilon\prime\prime_{m}=1$')
#pl.savefig('plots/131010_Hopkins_CG10_eps2_x_z.png', dpi = 300 )
pl.savefig('plots/131023_Hopkins_egn12_G_theta.pdf')
	#pl.axis([0.0,0.003,-0.1e-40,0.0])






