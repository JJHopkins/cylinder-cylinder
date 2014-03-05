#!/usr/bin/python
import numpy as np
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as pl

x_t_eV ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
x_x_eV ,y_x = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
x_y_eV ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
x_z_eV ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])

x_t =x_t_eV# /6.582e-16
x_x =x_x_eV# /6.582e-16
x_y =x_y_eV# /6.582e-16
x_z =x_z_eV# /6.582e-16

eiz_x = np.loadtxt('data/eiz_x_output.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output.txt')
eiz_z = np.loadtxt('data/eiz_z_output.txt') # parallel,axial
eiz_w = np.loadtxt('data/eiz_w_output.txt') # water as intervening medium
eiz_w[0] = eiz_w[1] #NOTE: there is a jump from first val down to second val
T = 300.0 
#kb = 8.6173e-5  # in eV/K
kb  = 1.3807e-23 # in J/K
c = 2.99e8

# NOTES:
# z_n_eV = (2*pi*kT/h_bar)n 
#	= (0.159 eV) / (6.5821e-16 eVs) 
#	= n*2.411e14 rad/s
# z_n_J = (2*pi*kT/h_bar)n 
#	= (1.3807e-23 J/K) / (1.0546e-34 Js))*n
#	= n*2.411e14 rad/s
#coeff = 0.159 # check me!
coeff = 2.411e14 # in rad/s
ns = np.arange(1.0,501.0)
z = ns * coeff

r_1 = 0.5e-9
r_2 = 0.5e-9
#ls = np.linspace(1e-9, 1000e-9, 10)
ls = np.linspace(1.0e-9, 7.0e-8, 10)
thetas = np.array([np.pi/5.0])
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
def ys(a,time,eizw,L, N):
	term0 = np.log( time    / (time*time+1.0)               )
	term1 =	np.log( time**4 * 2.0*(1. + 3.*a)*(1.+3.*a)     )
	term2 = np.log( time**2 * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a))
	term3 = np.log(          4.0*(1. + a)*(1.0 + a)         )
	term4 = (-2.0 * np.sqrt(eizw)* L * coeff * N * np.sqrt(time*time + 1.0)/(2.99e8))#*100))#NOTE: divided by 100 here
	#term4 = np.exp((-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8))
	#term4 = np.exp(-2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99)
#	term5 = 1.0 +    \
#		1e8*(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8) -    \
#		(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8) +    \
#		(1./2)*(1e8)**2*(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**2 +   \
#		(1./6)*(1e8)**3*(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**3
	#print 'ys term0', term0
	#print 'ys term1', term1
	#print 'ys term2', term2
	#print 'ys term3', term3
	#print 'ys term4', N,L,term4
#	rint 'ys term5', term5
#	print '----'
	#return term0 * term1 * term2 * term3 * term4 * term5
	return np.exp(term0 + term1 + term2 + term3 + term4)#* term5
#	return term0 * term1 * term2 * term3 * (1.0 + term4 +(1./2)*term4**2 +\
#		(1.0/6)*term4**3 + (1.0/24)*term4**4 +(1.0/120)*term4**5)

def y_2s(a,time,eizw, L, N): 
	term0 = np.log(time    / (time*time+1.0)                               ) 
	term1 = np.log((1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0))                     
	term2 = (-2.0*np.sqrt(eizw)*L*coeff*N*np.sqrt(time*time + 1.0)/(2.99e8))#*100)) #NOTE: divided by 100 here
	#term2 = np.exp((-2.0*L*coeff*N*np.sqrt(time*time + 1.0)/2.99e8))
#	term3 = np.log(1.0 +    \
#		1e8*(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8) -    \
#		(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8) +    \
#		(1./2)*(1e8)**2*(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**2 +   \
#		(1./6)*(1e8)**3*(2.0 * L * coeff * N * np.sqrt(time*time + 1.0)/2.99e8)**3
	#print 'y_2s term0', term0
	#print 'y_2s term1', term1
	#print 'y_2s term2', term2
#	print 'y_2s term3', term3
	#print '----'
	return np.exp(term0 + term1 + term2) #* term3


def As(eizz,eizw,L,N,Y): 
	term0 = (1.0/32)*kb*T
	term1 = ((eizz-eizw)/eizw)*((eizz-eizw)/eizw)
	term2 = eizw *eizw * (coeff*N)**4 * L**4 / ((2.99e8)**(4)) #NOTE: took out 1/c^4 for both A's
	term3 = Y
	#print 'As term2 with no e8= ', term2
	#print 'As term3 = ', term3
	return term0 * term1 * term2 * term3

def A_2s(eizz,eizw, L , N ,Y):		
	term0 = (1.0/32)*kb*T
	term1 = ((eizz-eizw)/eizw)*((eizz-eizw)/eizw)
	term2 = eizw *eizw * (coeff*N)**4 * L**4 * (2.99e8)**(-4)
	term3 = Y
	
	return term0 * term1 * term2 * term3

dt = 1000000000000000e-15
#ts = np.arange(1e-9,101e-9,dt)
ts = np.arange(1e-27,10000000000000000000e-15,dt)

A   = np.zeros(shape=(len(ns),len(ls)))
A_2 = np.zeros(shape=(len(ns),len(ls)))
aiz = []
sum_A = np.empty(len(ls))
sum_A_2 = np.empty(len(ls))
G = np.empty(len(thetas))
EL = np.zeros(len(ls))
aiz = Aiz(eiz_z,eiz_x, eiz_w) # of length = len(ns)
from pylab import pause
for j,n in enumerate(ns):
	print n
	for k,l in enumerate(ls):

		# Integrand:
		y_arg   = ys(aiz[j],ts,eiz_w[j],l,n)
		y_2_arg = y_2s(aiz[j],ts,eiz_w[j],l,n)
		
		#print 'ys = ',ys(aiz[j],ts,l,n)	
		#print 'y_2s = ',y_2s(aiz[j],ts,l,n)	
		#print '----'
		
		# Integral:
		y   = trapz(y_arg,ts,dt)
		y_2 = trapz(y_2_arg,ts,dt)
	#	print j,k	
	#	print 'int   y = ', y
	#	print 'int y_2 = ', y_2
	#	print 'As'  , As(eiz_z[j],l,n,y)
	#	print 'A_2s', A_2s(eiz_z[j],l,n,y_2)
	#	print '----'
		A[j,k]   = As(eiz_z[j],eiz_w[j],l,n,y)
		A_2[j,k] = A_2s(eiz_z[j],eiz_w[j],l,n,y_2) 		

sum_A = np.sum(A,axis=0)
#print 'sum_A = ', sum_A
sum_A_2 = np.sum(A_2,axis=0)
#print 'sum_A_2 = ', sum_A_2
#print '-----'

i = np.arange(len(ls))
EL[i] = 1./(ls[i]*ls[i]*ls[i]*ls[i])
G01 = -(EL[0]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[0] + sum_A_2[0]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G11 = -(EL[1]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[1] + sum_A_2[1]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G21 = -(EL[2]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[2] + sum_A_2[2]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G31 = -(EL[3]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[3] + sum_A_2[3]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G41 = -(EL[4]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[4] + sum_A_2[4]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G51 = -(EL[5]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[5] + sum_A_2[5]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G61 = -(EL[6]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[6] + sum_A_2[6]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G71 = -(EL[7]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[7] + sum_A_2[7]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G81 = -(EL[8]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[8] + sum_A_2[8]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))
G91 = -(EL[9]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[9] + sum_A_2[9]*np.cos(2.0*1.0/16*np.pi))/(2.0*np.sin(1.0/16*np.pi))

G02 = -(EL[0]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[0] + sum_A_2[0]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G12 = -(EL[1]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[1] + sum_A_2[1]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G22 = -(EL[2]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[2] + sum_A_2[2]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G32 = -(EL[3]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[3] + sum_A_2[3]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G42 = -(EL[4]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[4] + sum_A_2[4]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G52 = -(EL[5]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[5] + sum_A_2[5]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G62 = -(EL[6]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[6] + sum_A_2[6]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G72 = -(EL[7]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[7] + sum_A_2[7]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G82 = -(EL[8]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[8] + sum_A_2[8]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))
G92 = -(EL[9]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[9] + sum_A_2[9]*np.cos(2.0*1.0*np.pi/8))/(2.0*np.sin(1.0*np.pi/8))

G0 = -(EL[0]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[0] + sum_A_2[0]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G1 = -(EL[1]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[1] + sum_A_2[1]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G2 = -(EL[2]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[2] + sum_A_2[2]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G3 = -(EL[3]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[3] + sum_A_2[3]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G4 = -(EL[4]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[4] + sum_A_2[4]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G5 = -(EL[5]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[5] + sum_A_2[5]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G6 = -(EL[6]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[6] + sum_A_2[6]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G7 = -(EL[7]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[7] + sum_A_2[7]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G8 = -(EL[8]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[8] + sum_A_2[8]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))
G9 = -(EL[9]*np.pi*r_1*r_1*r_2*r_2)*(sum_A[9] + sum_A_2[9]*np.cos(2.0*3.0*np.pi/16))/(2.0*np.sin(3.0*np.pi/16))

print G01
print G11
print G21
print G31
print G41
print G51
print G61
print G71
print G81
print G91

print G02
print G12
print G22
print G32
print G42
print G52
print G62
print G72
print G82
print G92

print G0
print G1
print G2
print G3
print G4
print G5
print G6
print G7
print G8
print G9
Gs1 =[-G01,-G11,-G21,-G31,-G41,-G51,-G61,-G71,-G81,-G91]
Gs2 =[-G02,-G12,-G22,-G32,-G42,-G52,-G62,-G72,-G82,-G92]
Gs  =[-G0,-G1,-G2,-G3,-G4,-G5,-G6,-G7,-G8,-G9]
print Gs1 
print Gs2 
print Gs 
pl.figure()
pl.plot(ls,Gs1, label = r'$\theta = \frac{\pi}{16}$')
pl.plot(ls,Gs2, label = r'$\theta = \frac{\pi}{8}$')
pl.plot(ls,Gs , label = r'$\theta = \frac{3\pi}{16}$')
pl.legend()
pl.xlabel(r'$\mathrm{spearation}\,\it{l}\,\,\,\rm{[nm]}$', size = 20)
pl.ylabel(r'$\mathrm{G(\it{l})}\,\,\,\mathrm{[J]}$', size = 20)
pl.title(r'$\mathrm{G(\it{l},\theta= \pi/5)\,retarded,\,vacuum}$', size = 20)
pl.savefig('plots/130104_G_vs l_retarded.pdf' )
pl.show()
pl.figure()
pl.loglog(ls,Gs1, label = r'$\theta = \frac{\pi}{8}$')
pl.loglog(ls,Gs2, label = r'$\theta = \frac{2\pi}{8}$')
pl.loglog(ls,Gs , label = r'$\theta = \frac{3\pi}{8}$')
pl.legend()
pl.xlabel(r'$\mathrm{spearation}\,\it{l}\,\,\,\rm{[nm]}$', size = 20)
pl.ylabel(r'$\mathrm{G(\it{l})}\,\,\,\mathrm{[J]}$', size = 20)
pl.title(r'$\mathrm{G(\it{l},\theta= \pi/5)\,retarded,\,vacuum}$', size = 20)
pl.savefig('plots/130104_G_vs l_retarded.pdf' )
pl.show()
pl.close()


#for i in np.arange(0,10):
#	G = -(ls[i]**(-4)*np.pi*r_1*r_1*r_2*r_2)*(sum_A[i] +\
#	sum_A_2[i]*np.cos(2.0*np.pi/5))/(2.0*np.sin(np.pi/5))
#	print 'i,G,ls,A,A2 = ', i,G,ls[i],sum_A[i],sum_A_2[i]
###EL = ls**(-4)
###print 'EL =' , EL
###for i, el in enumerate(EL):
###	print 'el = ',el
###	G  = -(el*np.pi*r_1*r_1*r_2*r_2)*(sum_A[i] + sum_A_2[i]*np.cos(2.0*np.pi/5))/(2.0*np.sin(np.pi/5)) 
###	
###print G
#G=np.zeros(shape=(len(ls),len(thetas)))
#for t,theta in enumerate(thetas):
    #G[:,t] = -((L*np.pi*r_1*r_1*r_2*r_2)/(2.0*np.sin(theta)))\
#	*(sum_A + sum_A_2*np.cos(2.0*theta)) 
#G = -1.0*L*(sum_A + sum_A_2 )

lts = [':','-.','--','-']
markers = ['x', '^', 'd']
colors = ['b','g','r','c','y','m']

####pl.figure()
####pl.loglog(ls,-G)#,color = colors[j])
####pl.show()




