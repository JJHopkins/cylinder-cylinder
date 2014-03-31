#$ {\bf Free energy between two parallel cylinders (CG-10 in water). Nonretarded result, function of separation $\cal{l}$}
#$ Equation 31: $G(\ell,\theta) = - \frac{ (\pi R_1^{2})(\pi R_2^{2}) }{2 \pi~\ell^{4} \sin{\theta}} \left( {\cal A}^{(0)}(\ell) + {\cal A}^{(2)}(\ell) \cos 2\theta \right)$

#!/usr/bin/python
import numpy as np
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as pl
import pyreport
from matplotlib import axis as ax
# use pyreport -l file.py
from pylab import show
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from pylab import pause
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/par_NR_water/par_NR_water.pdf')

G_srw = np.loadtxt('g_prw.txt') 
eiz_x = np.loadtxt('data/eiz_z_output_eV.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output_eV.txt')
eiz_z = np.loadtxt('data/eiz_x_output_eV.txt') # parallel,axial
#eiz_w = 1.0 + np.zeros(len(eiz_z))
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # water as intervening medium
eiz_w[0] = eiz_w[1] #NOTE: there is a jump from first val down to second val

r_1 = 1.0e-9
r_2 = 1.0e-9
c = 2.99e8 # in m/s
Temp = 297 
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

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
def ys(a):
	term1 =	3.0 + 5.*(a + a)
	return term1
def y_2s(a): 
	term1 = (19.*a*a)                     
	return term1
def As(eizz,eizw,Y): 
	term1 = ((eizz-eizw)/eizw)*((eizz-eizw)/eizw)
	term2 = Y
	return term1 * term2
def A_2s(eizz,eizw, Y):		
	term1 = ((eizz-eizw)/eizw)*((eizz-eizw)/eizw)
	term2 = Y
	return term1 * term2

y = np.zeros(len(ns))#,len(ls)))
y_2 = np.zeros(len(ns))#,len(ls)))
A   = np.zeros(len(ns))
A_2 = np.zeros(len(ns))
sum_A = np.empty(len(ls))
sum_A_2 = np.empty(len(ls))
EL = np.zeros(len(ls))
G_l_t_dt = np.empty(len(ls))

aiz = []
aiz = Aiz(eiz_x,eiz_z, eiz_w) # of length = len(ns)

for j,n in enumerate(ns):
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
for k,length in enumerate(ls):
	EL[k] = 1./(length*length*length*length*length)
	G_l_t_dt[k] = (1.602e-19 / 4.11e-21) * (9./(32.*64.)) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A + sum_A_2)

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
#pl.savefig('plots/par_NR_water/eiz.pdf' )
#show()

pl.figure()
pl.loglog(ls, G_l_t_dt, label = r'$\mathrm{\mathcal{A}\,=\,}%3.2f$'%(1.0e21*(3.*kb*Temp/2)*(sum_A+sum_A_2)))
pl.xlabel(r'$\ell\,\,\mathrm{[m]}$', size = 24)
pl.ylabel(r'$\mathrm{-g(\ell,\theta=0;c\rightarrow\infty)\,\,[k_{B}T]}$', size = 20)
#pl.axis([1.9e-9, 8.1e-9,1e6,1e10])
#pl.title(r'$\mathrm{-g(\ell;c\rightarrow\infty)\,vs.\,separation:\,parallel,\,non-ret,\,water}$', size = 20)
pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/par_NR_water/140306_93w93_Nonret_skew_g_vs_l.pdf')
pl.show()

#ls_power = np.linspace(30.0e-9,1.0e-6,200)
ls_power = np.linspace(1.0e-9,1.0e-6,200)
ls_power_plus = 100.0e-9 + np.linspace(1.0e-9,1.0e-6,200)
l_5 = ls_power**(-5)
#l_5 = (30.0e-9 + ls_power)**(-5)
l_6 = (ls_power_plus)**(-6)
l_60 = (ls_power)**(-6)
#l_5 = (ls_power)**(-5)

pl.figure()
pl.loglog(ls_power,l_5/l_5[0],'b--', label = r'$\ell^{-5}$')
pl.loglog(ls_power_plus,l_6/(1.0e-2*l_60[0]),'m--', linewidth = 2.0, label = r'$\ell^{-6}$')

pl.loglog(ls, G_l_t_dt/G_l_t_dt[0],'k:', marker = 'o', mec = 'k', mew = 1, mfc = 'none', ms = 5, markevery = 5, label = r'$non-retarded$') 
pl.loglog(ls,    G_srw/   G_srw[0],'r:', marker = 's', mec = 'r', mew = 1, mfc = 'none', ms = 5, markevery = 5, label = r'$retarded$')
#pl.loglog(1.8e-9,0.9e-2,'b-',  label = r'$\theta = \pi/8$')
#pl.loglog(1.8e-9,0.9e-2,'g-',  label = r'$\theta = \pi/4$')
#pl.loglog(1.8e-9,0.9e-2,'r-',  label = r'$\theta = \pi/2$')
#pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
pl.xlabel(r'$\ell\,\,\mathrm{[m]}$', size = 24)
pl.ylabel(r'$-g(\ell,\theta=0)\,\,\mathrm{[k_{B}T]}$', size = 20)
pl.axis([1.0e-9, 1.0e-6,1.0e-16,1.0])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/par_NR_water/140306_93w93_Compare_par_normalized_g_vs_l.pdf')
pl.show()

