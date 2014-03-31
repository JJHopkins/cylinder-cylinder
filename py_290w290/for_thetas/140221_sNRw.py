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

G_srw = np.loadtxt('G_srw.txt') #perpendicular, radial

eiz_x = np.loadtxt('data/eiz_z_output_eV.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output_eV.txt')
eiz_z = np.loadtxt('data/eiz_x_output_eV.txt') # parallel,axial
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
ls = [3.0e-9,4.0e-9,10.0e-9,100.0e-9]
thetas = np.linspace((0.01)*np.pi,(1./2)*np.pi,25)
#ls = np.linspace(1.0e-9, 1.0e-6, 200)
#thetas = [np.pi/8,np.pi/4,np.pi/3,np.pi/2]
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
A0A2_theta = np.zeros(shape=(len(ls),len(thetas)))

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
sum_A_2 = np.sum(A_2,axis=0)
for k,length in enumerate(ls):
	for i, theta in enumerate(thetas):
		EL[k] = 1./(length*length*length*length)
		A0A2_theta[k,i] = sum_A + sum_A_2* np.cos(2.0*theta)
		G_l_t_dt[k,i] = (1.602e-19 / 4.11e-21) * (1./32) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A + sum_A_2* np.cos(2.0*theta) )/(2.0*np.sin(theta))# (1e21)*

#pl.figure()
#pl.loglog(ls,(kb*T/32)*sum_A,'b-', label = r'$\mathcal{A^{(0)}}(\ell)$')
#pl.loglog(ls,(kb*T/32)*sum_A_2,'g-', label = r'$\mathcal{A^{(2)}}(\ell)$')
##pl.loglog(ls,(kb*T/32)*A2_theta,':', label = r'$\mathcal{A^{(2)}}(\ell)cos(2\theta)$')
#pl.xlabel(r'$\mathrm{separation}\,\ell\,\,\,\rm{[m]}$', size = 20)
#pl.ylabel(r'$\mathrm{\mathcal{A^{(0)},\,\,A^{(2)}}}$', size = 20)
##pl.title(r'$\mathrm{Hamaker \, coeff.s \,:\,skewed,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'upper right')
##pl.axis([1e-9,1e-6,1e-24,1e-19])
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.savefig('plots/skew_NR_water/140306_290w290_skew_NR_A0_A2.pdf')
#show()

#ls1 = 1e9*ls[8]
#ls2 = 1e9*ls[12]
#ls3 = 1e9*ls[16]

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
#pl.semilogy(thetas, G_l_t_dt[0,:], 'b-', label = r'$\ell$ = %3.1f nm' %(1e9*ls[0]))
pl.semilogy(thetas, G_l_t_dt[1,:], 'b-', label = r'$\ell$ = %3.1f nm' %(1e9*ls[1]))
pl.semilogy(thetas, G_l_t_dt[2,:], 'g-', label = r'$\ell$ = %3.1f nm' %(1e9*ls[2]))
pl.semilogy(thetas, G_l_t_dt[3,:], 'r-', label = r'$\ell$ = %3.1f nm' %(1e9*ls[3]))

pl.semilogy(0, 1e-12, 'k-', label = r'$non-retarded$')
pl.semilogy(0, 1e-12, 'k:', label = r'$retarded$')

#pl.semilogy(thetas, G_srw[0,:], 'b--', linewidth = 2)#label = r'$\ell$ = %3.1f nm' %(1e9*ls[0]))
pl.semilogy(thetas, G_srw[1,:], 'b--', linewidth = 2)#label = r'$\ell$ = %3.1f nm' %(1e9*ls[1]))
pl.semilogy(thetas, G_srw[2,:], 'g--', linewidth = 2)#label = r'$\ell$ = %3.1f nm' %(1e9*ls[2]))
pl.semilogy(thetas, G_srw[3,:], 'r--', linewidth = 2)#label = r'$\ell$ = %3.1f nm' %(1e9*ls[3]))


#ax.semilogy(thetas, G_l_t_dt[ 2,:], label = r'$\ell$ = %1.2f nm' %ls4)
#ax.semilogy(thetas, G_l_t_dt[12,:], label = r'$\ell$ = %1.2f nm' %ls5)
##ax.semilogy(thetas, G_l_t_dt[22,:], label = r'$\ell$ = %1.2f nm' %ls6)
#ax.semilogy(thetas, G_l_t_dt[32,:], label = r'$\ell$ = %1.2f nm' %ls1)
##ax.semilogy(thetas, G_l_t_dt[42,:], label = r'$\ell$ = %1.2f nm' %ls2)
#ax.semilogy(thetas, G_l_t_dt[52,:], label = r'$\ell$ = %1.2f nm' %ls3)
##ax.semilogy(0,0,'', label = r'$G_\theta = cos(2\theta)/2sin(\theta)$')
pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
pl.axis([0,1.7,1e-11,1e3])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,retarded,\,water}$', size = 20)
pl.legend(loc = 'lower left')
#pl.savefig('plots/skew_ret_water/skew_ret_water_G_vs_theta.pdf')
#show()
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/skew_NR_water/140306_290w290_skew_compare_G_vs_theta.pdf')
show()

#pl.figure()
#pl.semilogy(thetas, G_l_t_dt[ 8,:], label = r'$\ell$ = %1.2f nm' %ls1)
#pl.semilogy(thetas, G_l_t_dt[12,:], label = r'$\ell$ = %1.2f nm' %ls2)
#pl.semilogy(thetas, G_l_t_dt[16,:], label = r'$\ell$ = %1.2f nm' %ls3)
#pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
#pl.ylabel(r'$-G(\ell,\theta;c\rightarrow\infty)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([(1./25)*np.pi,(3./4)*np.pi,105,135])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,non-ret,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.savefig('plots/skew_NR_water/Nonret_skew_G_vs_theta.pdf')
#show()
#
#Z = 1.0e21*(kb*T/32)*A0A2_theta[1,:]
#
#pl.figure()
##pl.loglog(ls, G_l_t_dt)#, label = labels[i])
#pl.loglog(ls, G_l_t_dt[:,1], label = r'$\theta = \pi/4,\,\mathcal{A} = %2.2f$'%Z[1])
#pl.loglog(ls, G_l_t_dt[:,2], label = r'$\theta = \pi/3,\,\mathcal{A} = %2.2f$'%Z[2])
#pl.loglog(ls, G_l_t_dt[:,3], label = r'$\theta = \pi/2,\,\mathcal{A} = %2.2f$'%Z[3])
#pl.loglog(ls[0],G_l_t_dt[0,0],' ', label = r'$\mathcal{A}^{(0)} = %3.2f$'%(1.0e21*(kb*T/32)*sum_A))
#pl.loglog(ls[0],G_l_t_dt[0,0],' ', label = r'$\mathcal{A}^{(2)} = %3.2f$'%(1.0e21*(kb*T/32)*sum_A_2))
##pl.loglog(ls, (kb*T/32)*A0A2_theta[:,1],'g ', label = r'$\mathcal{A}(\ell,\theta = \pi/4) = %2.2f$'%Z[1])#G_l_t_dt[:,3], label = r'$\theta = \pi/2$')
#
#pl.xlabel(r'$\ell\,\,\mathrm{[m]}$', size = 24)
#pl.ylabel(r'$-G(\ell,\theta;c\rightarrow\infty)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([1.0e-9, 1.0e-6,1e-12,1e2])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.savefig('plots/skew_NR_water/140306_290w290_Nonret_skew_G_vs_l.pdf')
#show()
#
##ls_power = np.linspace(30.0e-9,1.0e-6,200)
#ls_power = np.linspace(1.0e-9,1.0e-6,200)
#ls_power_plus = 100.0e-9 + np.linspace(1.0e-9,1.0e-6,200)
#l_4 = ls_power**(-4)
##l_5 = (30.0e-9 + ls_power)**(-5)
#l_5 = (ls_power_plus)**(-5)
#l_50 = (ls_power)**(-5)
##l_5 = (ls_power)**(-5)
#
#pl.figure()
#pl.loglog(ls_power,l_4/l_4[0],'b--', label = r'$\ell^{-4}$')
#pl.loglog(ls_power_plus,l_5/(1.0e-2*l_50[0]),'m--', linewidth = 2.0, label = r'$\ell^{-5}$')
#pl.loglog(ls, G_l_t_dt[:,3]/G_l_t_dt[0,3],'k:', marker = 'o', mec = 'k', mew = 1, mfc = 'none', ms = 5, markevery = 5, label = r'$non-retarded$') 
#pl.loglog(ls,    G_srw[:,3]/   G_srw[0,3],'r:', marker = 's', mec = 'r', mew = 1, mfc = 'none', ms = 5, markevery = 5, label = r'$retarded$')
##pl.loglog(1.8e-9,0.9e-2,'b-',  label = r'$\theta = \pi/8$')
##pl.loglog(1.8e-9,0.9e-2,'g-',  label = r'$\theta = \pi/4$')
##pl.loglog(1.8e-9,0.9e-2,'r-',  label = r'$\theta = \pi/2$')
##pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
#pl.xlabel(r'$\ell\,\,\mathrm{[m]}$', size = 24)
#pl.ylabel(r'$-G(\ell,\theta=\pi/2)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.ylabel(r'$-G(\ell,\theta=\frac{\pi}{2})\,\,\mathrm{[k_{B}T]}$', size = 20)
#pl.axis([1.0e-9, 1.0e-6,1.0e-13,1.0])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
##pl.savefig('plots/Compare_full_NR_G_vs_l_normalized.pdf')
#pl.savefig('plots/skew_NR_water/140306_290w290_Compare_normalized_G_vs_l_perp.pdf')
#show()
#
