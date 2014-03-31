#$ {\bf Free energy between two skewed cylinders (CG-10 in water). Nonretarded result, function of separation $\ell$ and angle $\theta$} \\
#$ Equation 18: $G(\ell,\theta; c \longrightarrow \infty) &=&  - \frac{k_BT}{64 \pi} \frac{\pi^{2} R_1^{2} R_2^{2}}{\ell^{4} \sin{\theta}} {\sum_{n=0}^{\infty}}' \Delta_{1,\parallel} \Delta_{2,\parallel} \frac{3}{8}\left[ 2 (1+3a_1) (1+3a_2) + (1-a_1) (1-a_2)  \cos 2\theta \right].$

#!/usr/bin/python
import numpy as np
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as pl
from matplotlib import axis as ax
# use pyreport -l file.py
from pylab import show
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from pylab import pause
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/skew_NR_water/skew_NR_water.pdf')

Gsrw = np.loadtxt('G_srw.txt') # retarded skew output for len(ls,theta) = 5,5 

eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output_eV.txt')
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # parallel,axial
#eiz_w = 1.0 + np.zeros(len(eiz_z))
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # water as intervening medium
eiz_w[0] = eiz_w[1] #NOTE: there is a jump from first val down to second val

r_1 = 1.0e-9
r_2 = 1.0e-9
c = 2.99e8 # in m/s
T = 297 
kb  = 1.3807e-23 # in J/K
coeff = 2.411e14 # in rad/s
# NOTES:
# h_bar = 1. #1.0546e-34 #in Js
#kb = 8.6173e-5  # in eV/K
# at RT, 1 kT = 4.11e-21 J
# 1 eV = 1.602e-19 J = 0.016 zJ
# h_bar_eV = 6.5821e-16 eVs
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
thetas = np.linspace((0.01)*np.pi,(1./2)*np.pi,25)
dt = 1.0
ts = np.arange(1.0,10000.,dt)

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
def ys(a):
	term1 =	(2.0*(1. + 3.*a)*(1.+3.*a))
	return term1
def y_2s(a): 
	term1 = ((1.- a)*(1.- a))                     
	return term1
def As(eizz,eizw,Y): 
	term1 = 3./8*((eizz-eizw)/eizw)*((eizz-eizw)/eizw)
	term2 = Y
	return term1 * term2
def A_2s(eizz,eizw,Y):		
	term1 = 3./8*((eizz-eizw)/eizw)*((eizz-eizw)/eizw)
	term2 = Y
	return term1 * term2

y = np.zeros(len(ns))#,len(ls)))
y_2 = np.zeros(len(ns))#,len(ls)))
A   = np.zeros(len(ns))#,len(ls)))
A_2 = np.zeros(len(ns))#,len(ls)))
EL = np.zeros(len(ls))
G_l_t_dt = np.zeros(shape=(len(ls),len(thetas)))

aiz = []
aiz = Aiz(eiz_x,eiz_z, eiz_w) # of length = len(ns)

sum_A = np.empty(len(ls))
sum_A_2 = np.empty(len(ls))

for j,n in enumerate(ns):
	# Integral:
	y[j]   = ys(aiz[j])
	y_2[j] = y_2s(aiz[j])
	#print 'dt Integral   y = ',i,k,j, y
	#print 'dt Integral y_2 = ',i,k,j, y_2
	#print '----'
	#print 'N terms for A0 = '  , As(eiz_z[j],eiz_w[j],length,n,y)
	#print 'N terms for A2 = ', A_2s(eiz_z[j],eiz_w[j],length,n,y_2)
	#print '----'
	A[j]   = As(eiz_z[j],eiz_w[j],y[j])
	A_2[j] = A_2s(eiz_z[j],eiz_w[j],y_2[j])# * np.cos(2.0*theta)  		
	A[0] = (1./2)*A[0]  		
	A_2[0] = (1./2)*A_2[0]  		
sum_A = np.sum(A,axis=0)
#print 'sum of A0 = ', j,sum_A
sum_A_2 = np.sum(A_2,axis=0)
#print 'sum of A2 = ', j,sum_A_2
#print '----'
#print 'shape sum_A_2 = ', np.shape(sum_A_2)
#sys.exit()
for k,length in enumerate(ls):
	for i, theta in enumerate(thetas):
		EL[k] = 1./(length*length*length*length)
		G_l_t_dt[k,i] = (1.602e-19 / 4.11e-21) * (1./32) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A + sum_A_2* np.cos(2.0*theta) )/(2.0*np.sin(theta))# (1e21)*

ls1 = 1e9*ls[ 2]
ls2 = 1e9*ls[12]
ls3 = 1e9*ls[22]
ls4 = 1e9*ls[32]
ls5 = 1e9*ls[42]
ls6 = 1e9*ls[52]

#pl.figure()
#pl.loglog(ls, G_l_t_dt,'-', ls, Gsrw, ':')#, label = labels[i])
#pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
#pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([1.5e-9, 6.5e-8,100,145])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
##pl.legend(loc = 'best')
##pl.savefig('plots/skew_NR_water/G_vs_l.pdf')
#show()
ls2 = np.linspace(1e-8,1e-7,25)
l_4 =ls2**(-4)# (1.602e-19 / 4.11e-21) * (1./32) *ls**(-4)
l_5 =ls2**(-5)# (1.602e-19 / 4.11e-21) * (1./32) *ls**(-5)

pl.figure()
pl.semilogy(ls, G_l_t_dt[:, 6],'b-', ls, Gsrw[:, 6], 'b:')#, label = labels[i])
pl.semilogy(ls, G_l_t_dt[:,12],'g-', ls, Gsrw[:,12], 'g:')#, label = labels[i])
pl.semilogy(ls, G_l_t_dt[:,16],'r-', ls, Gsrw[:,16], 'r:')#, label = labels[i])
pl.semilogy(1.8e-9,0.9e-2,'b-',  label = r'$\theta = \pi/8$')
pl.semilogy(1.8e-9,0.9e-2,'g-',  label = r'$\theta = \pi/4$')
pl.semilogy(1.8e-9,0.9e-2,'r-',  label = r'$\theta = \pi/3$')
pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
pl.axis([1.0e-9, 1.0e-6,1e-12,1e2])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
pl.legend(loc = 'best')
#pl.savefig('140220_cyl-cyl/Compare_full_NR_G_vs_l.pdf')
show()

pl.figure()
pl.loglog(ls, G_l_t_dt[:, 6],'b-', ls, Gsrw[:, 6], 'b:')#, label = labels[i])
pl.loglog(ls, G_l_t_dt[:,12],'g-', ls, Gsrw[:,12], 'g:')#, label = labels[i])
pl.loglog(ls, G_l_t_dt[:,16],'r-', ls, Gsrw[:,16], 'r:')#, label = labels[i])
pl.loglog(1.8e-9,0.9e-2,'b-',  label = r'$\theta = \pi/8$')
pl.loglog(1.8e-9,0.9e-2,'g-',  label = r'$\theta = \pi/4$')
pl.loglog(1.8e-9,0.9e-2,'r-',  label = r'$\theta = \pi/3$')
pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
pl.axis([1.0e-9, 1.0e-6,1e-12,1e3])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
pl.legend(loc = 'best')
pl.savefig('140220_cyl-cyl/Compare_full_NR_G_vs_l.pdf')
show()

pl.figure()
#pl.loglog(ls2,l_4/l_4[0],'k-' )#G_l_t_dt[:, 0],'b-', ls, Gsrw[:, 0], ':')#, label = labels[i])
pl.loglog(ls2,l_4,'k-' )#G_l_t_dt[:, 0],'b-', ls, Gsrw[:, 0], ':')#, label = labels[i])
#pl.loglog(ls2,l_5/l_5[0],'k:' )#G_l_t_dt[:, 5],'b-', ls, Gsrw[:, 5], 'b:')#, label = labels[i])
pl.loglog(ls2,l_5,'k:' )#G_l_t_dt[:, 5],'b-', ls, Gsrw[:, 5], 'b:')#, label = labels[i])
pl.loglog(ls, G_l_t_dt[:, 6]/G_l_t_dt[0, 6],'b-', ls, Gsrw[:, 6]/Gsrw[0, 6], 'b--')#, label = labels[i])
pl.loglog(ls, G_l_t_dt[:,12]/G_l_t_dt[0,12],'g-', ls, Gsrw[:,12]/Gsrw[0,12], 'g--')#, label = labels[i])
pl.loglog(ls, G_l_t_dt[:,16]/G_l_t_dt[0,16],'r-', ls, Gsrw[:,16]/Gsrw[0,16], 'r--')#, label = labels[i])
pl.loglog(1.8e-9,0.9e-2,'b-',  label = r'$\theta = \pi/8$')
pl.loglog(1.8e-9,0.9e-2,'g-',  label = r'$\theta = \pi/4$')
pl.loglog(1.8e-9,0.9e-2,'r-',  label = r'$\theta = \pi/3$')
pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
pl.axis([1.0e-9, 1.0e-6,1e-12,1e3])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
pl.legend(loc = 'best')
pl.savefig('140220_cyl-cyl/Compare_full_NR_G_vs_l_normalized.pdf')
show()

pl.figure()
pl.semilogy(thetas, G_l_t_dt[ 2,:],'b-',thetas , Gsrw[ 2,:], 'b:')#, label = labels[i]), label =r'$\ell$ = %1.2f nm' %ls1
pl.semilogy(thetas, G_l_t_dt[12,:],'g-',thetas , Gsrw[12,:], 'g:')#, label = labels[i]), label =r'$\ell$ = %1.2f nm' %ls2
pl.semilogy(thetas, G_l_t_dt[22,:],'r-',thetas , Gsrw[22,:], 'r:')#, label = labels[i]), label =r'$\ell$ = %1.2f nm' %ls3
pl.semilogy(thetas, G_l_t_dt[32,:],'c-',thetas , Gsrw[32,:], 'c:')#, label = labels[i]), label =r'$\ell$ = %1.2f nm' %ls4
pl.semilogy(thetas, G_l_t_dt[42,:],'m-',thetas , Gsrw[42,:], 'm:')#, label = labels[i]), label =r'$\ell$ = %1.2f nm' %ls5
pl.xlabel(r'$Angle\,\,\theta\,\,\mathrm{[radians]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
#pl.axis([1.5e-9, 6.5e-8,100,145])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.savefig('plots/skew_NR_water/G_vs_l.pdf')
show()



#G_l_t_dt[G_l_t_dt>300]= np.nan #NOTE: remove me later
#G_l_t_dt[G_l_t_dt<200e-25]= np.nan #NOTE: remove me later
