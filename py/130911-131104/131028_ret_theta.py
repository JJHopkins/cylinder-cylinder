#!/usr/bin/python

import numpy as np
from numpy import log,exp
#from scipy import *
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import axis as ax
import matplotlib               
from pylab import *
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

x_t_eV ,y_t  = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x_eV ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y_eV ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z_eV ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

x_t =x_t_eV /6.582e-16
x_x =x_x_eV /6.582e-16
x_y =x_y_eV /6.582e-16
x_z =x_z_eV /6.582e-16

eiz_x = np.loadtxt('data/eiz_x_output.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') # parallel,axial

coeff = 2.41 # in rad/s
ns = arange(1.0,500.0)
T = 300.0 # r.t. in K
kb_J =1.3806488e-23# in zJ/K # 1.3806488e-23 # in J/K
hbar =6.582e-16 # in inverse sec #6.625e-34 # in J/s
z = ns * coeff
r_1 = 0.5e-9
r_2 = 0.5e-9
#r_1 = 0.5
#r_2 = 0.5
c = 2.99e8

#t = np.linspace(1e-9,1e-7,10)
#ls = [1e-9 , 50e-9 , 100e-9]
#ls = linspace(0.1e-9, 10e-9, 100)
ls = linspace(1.1e-9, 1000e-9, 10)
#ls = [4.0]
#thetas = linspace(0.0,np.pi/2.0,100)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
#thetas = linspace(np.pi/5.0,np.pi/2.0,100)#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
thetas = [np.pi/5.0,np.pi/4.0,np.pi/3.0]#np.pi, 2*np.pi,6) #/4, np.pi/3, np.pi][np.pi/2.0,2.0*np.pi]#
eiz_m = 1.0

#height= np.zeros(shape=(len(ls),len(thetas)))
lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']
#
aiz = []#np.zeros(len(eiz_z))
#A = np.empty(len(ns))
#A_2 = np.empty(len(ns))
#y = np.empty(len(ns))
#y_2 = np.empty(len(ns))
#g = np.empty(len(ns))
#g_int = np.empty(len(ns))
#A0 = np.empty(len(ns))
#w = np.empty(len(ls))
sum_A = np.empty(len(ls))
sum_A_2 = np.empty(len(ls))
#sum_A_thetas = np.empty(len(thetas))
#g_t = np.empty(len(thetas))
#A_t = np.empty(len(thetas))
#G = np.empty(len(ls))
#added = np.empty(len(ls))
#added_thetas = np.empty(len(thetas))
#
def Aiz(perp, par):
	return 2.0*(perp-1.0)/((perp+1.0)*(par-1.0))

def ys(a,time,L, N):
	return (\
		(((time**5/(time*time+1.0)) * 2.0*(1. + 3.*a)*(1.+3.*a)     ) + \
		((time**3/(time*time+1.0)) * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a)) + \
		((time   /(time*time+1.0)) * 4.0*(1. + a)*(1.0 + a))) * \
		np.exp((-2.0)*L*(1./(2.99e8))*coeff*N*(time*time + 1.0)**(1./2)))
def y_2s(a,time, L, N): 
	return (\
		((time / (time * time + 1.0)) * (1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0)) * \
		np.exp((-2.0)*L*(1./(2.99e8))*coeff*N*(time*time  + 1.0)**(1./2)))
#def ys(a,time,L, N):
#	return (\
#		log((time**5/(time*time+1.0)) * 2.0*(1. + 3.*a)*(1.+3.*a)     ) + \
#		log((time**3/(time*time+1.0)) * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a)) + \
#		log((time   /(time*time+1.0)) * 4.0*(1. + a)*(1.0 + a)) + \
#		3.0*(-2.0)*L*(1./(2.99e8))*coeff*N*(time*time + 1.0)**(1./2))
#def y_2s(a,time, L, N): 
#	return (\
#		log((time / (time * time + 1.0)) * (1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0)) + \
#		3.0*(-2.0)*L*(1./(2.99e8))*coeff*N*(time*time  + 1.0)**(1./2))
def As(eizz,L,N,Y): 
	return (1.0/32)*kb_J*T*((eizz-1.0)/1.0)*((eizz-1.0)/1.0)*\
		1.0*((coeff*N)*(coeff*N)*L*L/((2.99e8)*(2.99e8)))*((coeff*N)*(coeff*N)*L*L/((2.99e8)*(2.99e8)))*\
		Y
		#np.exp(Y)
def A_2s(eizz, L , N ,Y):		
	return (1.0/32)*kb_J*T*((eizz-1.0)/1.0)*((eizz-1.0)/1.0)*\
		1.0*((coeff*N)*(coeff*N)*L*L/((2.99e8)*(2.99e8)))*((coeff*N)*(coeff*N)*L*L/((2.99e8)*(2.99e8)))*\
		Y
		#np.exp(Y)
#for i,ezz in enumerate(eiz_z):
#	for m,ezx in enumerate(eiz_x):
#		aiiz = Aiz(ezz,ezx)
#	aiz.append(aiiz)
aiz = np.zeros(len(eiz_z))
for i,ezz in enumerate(eiz_z):
	aiz = [Aiz(ezx,ezz) for ezx,ezz in zip(eiz_x,eiz_z)]
#aiz.append(aiiz)

ts = linspace(1e-9,1000e-1)
dt =10e-9
pl.figure()
for j,theta in enumerate(thetas):
	for p, l in enumerate(ls):
		print theta,l
		y_arg   = [ys(az,ts,l,n) for az,n in zip(aiz,ns)]
		y_2_arg = [y_2s(az,ts,l,n)for az,n in zip(aiz,ns)]
		y = trapz(y_arg,ts,dt)#axis = 0)
		y_2 = trapz(y_2_arg,ts,dt)#,axis =0)
		A   = [As(ezz,l,n,yy)   for ezz,yy in zip(eiz_z,y)]
		#print l,A
		A_2 = [A_2s(ezz,l,n,yy_2) for ezz,yy_2 in zip(eiz_z,y_2)]		
		sum_A[p] = sum(A)
		sum_A_2[p] = sum(A_2)
		#print sum_A 
		G = -((pi*r_1*r_1*r_2*r_2)/(l*l*l*l*2.0*sin(theta)))*(sum_A + sum_A_2*cos(2.0*theta)) 
	pl.loglog(ls,-G,color = colors[j])
pl.show()
##for j,n in enumerate(ns):
#	ts = linspace(1e-9,101e-9)
#	dt =10e-9
#	for k,ezz in enumerate(eiz_z):
#		for p, l in enumerate(ls):
#			for q,az in enumerate(aiz):
#				y_arg   = [ys(az,ts,l,n)]
#				y_2_arg = [y_2s(az,ts,l,n)]
#				y = trapz(y_arg,ts,dt)#axis = 0)
#				y_2 = trapz(y_2_arg,ts,dt)#,axis =0)
#			A   = [As(ezz,l,n,yy)   for yy in y]
#			A_2 = [A_2s(ezz,l,n,yy_2) for yy_2 in y_2]		
#
#			print'n_aiz =',q	
#	print 'n=',n
##	#				for p,theta in enumerate(thetas):
#		
#


#
#pl.figure()
#for i,phi_b in enumerate(phi_bs):
#	phis = np.linspace(0.0, 1.0-phi_b-phi_water_min, 100)
##	pl.subplot(len(phi_bs),1,i+1)
##	for df in dfs:
#	for j,df in enumerate(dfs):
#		try: ps = [P(phi,phi_b,df) for phi in phis]
#	       	except: continue
#		resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
#		max_resid = max(resids)
#		print 'Largest residual for df=%f: %e' % (df, max_resid)
#		if j==3:
#			labels=r'$\Phi_{PEG3500} = %.2f$' % phi_b
#		else:
#			labels=None
#		pps = [phi*P(phi,phi_b,df) for phi in phis]
#		print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
#		pl.plot(phis,pps,color=colors[j],linestyle = lts[i],label=labels)#color = colors[j],
##
##		#pl.text(0.9,1.95,r'$\Delta F = %.2f$' % j(0), color='r')
##		#pl.text(0.9,1.85,r'$\Delta F = %.2f$' % j(1), color='g')
##		#pl.text(0.9,1.75,r'$\Delta F = %.2f$' % j(2), color='b')
##		#pl.text(0.9,1.65,r'$\Delta F = %.2f$' % j(3), color= 'k')
##
##
##	pl.axis([0, 1.05, -0.1, 1.05])
##
###font = {'size': 22}
###matplotlib.rc('font',**font)
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
#
#
#
#
#
#
#
#
#
#
#
##			for i in range(len(ts)):
#			#A = ((eiz_z-1)/1)*((eiz_z-1)/1)*(2*((1 +3 * aiz)*(1+3 * aiz))+(1 - aiz)(1 - aiz)*cos(thetas[p]))
#			#sum_A[k] = sum(A)
#			#sum_A_thetas[p] = sum(A[j])
#
#		#	A = np.nan_to_num(A)
#			#A = [(3.0/8.0)*((eiz_zs-1.0)/1.0)*((eiz_zs-1.0)/1.0)*(2.0*((1.+3.*aizs)*(1.+3.*aizs))+(1.-aizs)*(1.0 - aizs)*cos(thetas[p])) for eiz_zs,aizs in zip(eiz_z,aiz)]
#			#g[j] = 2.0*((1.0+3.0*aiz[j])*(1.0+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p])
#
#			##g[j] = trapz(t/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c)*\
#			##(2.0*((1.+3.*aiz[j])*(1.+3.*aiz[j])*t*t*t*t+2.*(1.+2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])*t*t+2.0*(1.0+aiz[j])*(1.0 + aiz[j]))),t)
#			y[j] = trapz(\
#					log((t**5/(t*t+1.0)) * 2.0*(1. + 3.*aiz[j])*(1.+3.*aiz[j])               ) + \
#					log((t**3/(t*t+1.0)) * 4.0*(1. + 2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])) + \
#					log((t   /(t*t+1.0)) * 4.0*(1. + aiz[j])*(1.0 + aiz[j])) + \
#					3.0*(-2.0)*ls[k]*(1./c)*coeff*ns[j]*(t*t + 1.0)**(1./2), t)
#			 
#			y_2[j] = trapz(\
#					log((t/(t*t+1.0)) * (1.- aiz[j])*(1.- aiz[j])*(t*t + 2.0)*(t*t + 2.0)) + \
#					3.0*(-2.0)*ls[k]*(1./c)*coeff*ns[j]*(t*t + 1.0)**(1./2), t)
#			#z[j] = 3.0*2.0*ls[k]*(1./c)*coeff*ns[j] 
#
#			#g[j] =trapz(np.exp( log(t/(t*t+1.0)) + (-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j]*ls[k]/c) +\
#			#log(2.0*((1.+3.*aiz[j])*(1.+3.*aiz[j])*t*t*t*t+2.*(1.+2.0*aiz[j]+2.0*aiz[j]+3.0*aiz[j]*aiz[j])*t*t+2.0*(1.0+aiz[j])*(1.0 + aiz[j])))),t)
#			
#			#g_int[j] = trapz(g,t)#/(t*t+1.0)* np.exp(-2.0*1.0*(t*t+1.0)**(1./2)*coeff*ns[j])*\
#
#			#A[j] = (3.0/8.0)*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*(2.0*((1.+3.0*aiz[j])*(1.+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p]))
#
#			#A[j] = kb_J*T/32.0*((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*\
#			#1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
#			#exp(y[j])
#			A[j] = ((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*\
#				1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
#				exp(y[j])
#	
#			A_2[j] = ((eiz_z[j]-1.0)/1.0)*((eiz_z[j]-1.0)/1.0)*\
#				1.0*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*((coeff*ns[j])*(coeff*ns[j])*ls[k]*ls[k]/(c*c))*\
#				exp(y_2[j])
#			#A_t[k] = sum(A[j])#+= 2.0*((1.0+3.0*aiz[j])*(1.0+3.0*aiz[j]))+(1.0-aiz[j])*(1.0 - aiz[j])*cos(2.0*thetas[p])
#		pl.plot(ns,y)
#		sum_A[k] = sum(A)
#		sum_A_2[k] = sum(A_2)
#		w[k] = (r_1*r_1*r_2*r_2)/(ls[k]*ls[k]*ls[k]*ls[k])
#	#G[k] = -(np.pi*np.pi*r_1*r_1*r_2*r_2)/(np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p]))*sum_A[k]
#		zip(sum_A,sum_A_2)
#		G[k] = -(pi/sin(thetas[p]))*(sum_A[k] + sum_A_2[k] * cos(2.0*thetas[p])) 
#		q = 1.0/ls**4
#		Q = r_1**4
#		qs= q*Q
#		qs[0] = qs[2]
#		qs[1] = qs[2]
#		qs[49] = qs[48]
#		QS = qs*G
#	
#	pl.plot(ls,QS, color = colors[p])#, linestyle = lts[k])	
#pl.title(r'G vs l', size = 16)
#pl.figure()
#pl.plot(ns,y)
#pl.title(r'y vs n', size = 16)
#pl.figure()
#plot(ns,exp(y))			#print sum_A
#pl.title(r'exp(y) vs n', size = 16)
#			#sum_A_thetas[p] = sum(A[j])
#		#G = -(1.6e-19*(kb*T)*(np.pi*np.pi*r_1*r_1*r_2*r_2)/(64.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A
#		#G[k] = -((kb_J*T)*(np.pi*np.pi*r_1*r_1*r_2*r_2)/(64.0*np.pi*ls[k]*ls[k]*ls[k]*ls[k]*sin(thetas[p])))*Sum_A
##pl.figure()
#	#pl.plot(ns,g, color = colors[p])#, linestyle = lts[k])	
##pl.plot(ls,G, color = colors[p])#, linestyle = lts[k])	
#	#pl.semilogy(ls,-G, color = colors[p])#, linestyle = lts[k])	
##	pl.figure()
#
##pl.figure()
##pl.plot(ls, 26.0)#, color = colors[p])#, linestyle = lts[k])	
##pl.xlabel(r'$\theta\,\,\,[radians]$', size = 22)
##pl.ylabel(r'$ g(\theta)\,-\,26\,\rm{zJ/nm}$', size = 22)
##pl.title(r'g(l = 4 nm, $\theta$) Non retarded skewed cylinders', size = 16)
##pl.savefig('plots/131027_Hopkins_ret_g_theta.pdf')
##
##pl.figure()
##pl.plot(ls,A_t)#, color = colors[p])#, linestyle = lts[k])	
##pl.xlabel(r'$\theta\,\,\,[radians]$', size = 22)
##pl.ylabel(r'$ A(\theta)/A(0)$', size = 22)
##pl.title(r'$\Sigma \Delta_{1} \Delta_{2} \times \rm{g(l = 4 nm, \theta)}$ for non retarded skewed cylinders', size = 16)
##pl.savefig('plots/131027_Hopkins_ret_sum_theta.pdf')
#
#
#	#pl.text(0.1,-0.15e-8,r'$\theta = %.2f \, radians$' % thetas[0], color='r')
#	#pl.text(0.1,-0.30e-8,r'$\theta = %.2f \, radians$' % thetas[1], color='g')
#	#pl.text(0.1,-0.45e-8,r'$\theta = %.2f \, radians$' % thetas[2], color='b')
##pl.legend()
##pl.savefig('plots/131010_Hopkins_CG10_eps2_x_z.png', dpi = 300 )
#pl.show()
#pl.close()
#
#
#
#
#
#
#
#
