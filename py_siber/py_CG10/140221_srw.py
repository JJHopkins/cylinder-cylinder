#$ {\bf Free energy between two skewed cylinders (CG-10 in water). Full retarded result, function of separation $\ell$ and angle $\theta$} \\
#$ Equation 12: $G(\ell,\theta) = - \frac{ (\pi R_1^{2})(\pi R_2^{2}) }{2 \pi~\ell^{4} \sin{\theta}} \left( {\cal A}^{(0)}(\ell) + {\cal A}^{(2)}(\ell) \cos 2\theta \right)$ \\
#$ $G(\ell,\theta) = - \frac{k_BT}{64 \pi} \frac{ \pi^2 R_1^{2} R_2^{2} }{\ell^{4} \sin{\theta}} {\sum_{n=0}^{\infty}}' \Delta_{1,\parallel} \Delta_{2,\parallel} ~p_n^{4} ~\int_0^{\infty} t dt ~\frac{e^{- 2 p_n \sqrt{t^{2} + 1}}}{(t^{2} + 1)} \tilde g(t, a_1(i \omega_n), a_2(i \omega_n), \theta),$ \\
#$ with $\tilde g(t, a_1, a_2, \theta) &=& 2 \left[ (1+3a_1)(1+3a_2) t^{4} + 2 (1+2a_1+2a_2+3a_1a_2) t^{2}  + 2(1+a_1)(1+a_2)\right] + \nonumber \\
#$ & & ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ + (1-a_1)(1-a_2)(t^{2} + 2)^2 \cos 2\theta.$ \\

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
pp = PdfPages('plots/skew_ret_water/skew_ret_water.pdf')

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
#thetas = np.linspace((0.01)*np.pi,(1./2)*np.pi,25)
thetas = [np.pi/8,np.pi/4,np.pi/3,np.pi/2]
dt = 1.0
ts = np.arange(1.0,10000.,dt)

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
def ys(a,time,eizw,L, N):
	term0 = ( time    / (time*time+1.0)               )
	term1 =	( time**4 * 2.0*(1. + 3.*a)*(1.+3.*a)     )
	term2 = ( time**2 * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a))
	term3 = (          4.0*(1. + a)*(1.0 + a)         )
	term4 = (-2.0 * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	#print 'ys term0', term0
	#print 'ys term1', term1
	#print 'ys term2', term2
	#print 'ys term3', term3
	#print 'ys term4', term4
	#print '----'
	return (term0) * np.exp(term4)*( (term1) + (term2) + (term3))#* term5
def y_2s(a,time,eizw, L, N): 
	term0 = (time    / (time*time+1.0)                               ) 
	term1 = ((1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0))                     
	term2 = (-2.0 * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	#print 'y_2s term0', term0
	#print 'y_2s term1', term1
	#print 'y_2s term2', term2
	#print '----'
	return term0 * term1* np.exp(term2) #* term3
def As(eizz,eizw,L,N,Y): 
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * eizw *eizw * (coeff*N)**4 * L**4 / (c**4))
	#term3 = Y
	#print 'As term1 = ', term1
	#print 'As term2 = ', term2
	##print 'As term3 = ', term3
	#print '----'
	return  term1 *  term2# * term3
def A_2s(eizz,eizw, L , N ,Y):		
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * eizw *eizw * (coeff*N)**4 * L**4 / (c**4))
	#term3 = Y
	#print 'A_2s term1 = ', term1
	#print 'A_2s term2 = ', term2
	##print 'A_2s term3 = ', term3
	#print '----'
	return (term1 * term2)# * term3


y = np.zeros(shape=(len(ns),len(ls)))
y_2 = np.zeros(shape=(len(ns),len(ls)))
A   = np.zeros(shape=(len(ns),len(ls)))
A_2 = np.zeros(shape=(len(ns),len(ls)))
EL = np.zeros(len(ls))
G_l_t_dt = np.zeros(shape=(len(ls),len(thetas)))

aiz = []
aiz = Aiz(eiz_x,eiz_z, eiz_w) # of length = len(ns)

for k,length in enumerate(ls):
	sum_A = np.empty(len(ls))
	sum_A_2 = np.empty(len(ls))
	for j,n in enumerate(ns):
		# Integral:
		y[j,k]   = trapz(ys(aiz[j],ts,eiz_w[j],length,n),ts,dt)
		y_2[j,k] = trapz(y_2s(aiz[j],ts,eiz_w[j],length,n),ts,dt)
		#print 'dt Integral   y = ',i,k,j, y
		#print 'dt Integral y_2 = ',i,k,j, y_2
		#print '----'
		#print 'N terms for A0 = '  , As(eiz_z[j],eiz_w[j],length,n,y)
		#print 'N terms for A2 = ', A_2s(eiz_z[j],eiz_w[j],length,n,y_2)
		#print '----'
		A[j,k]   = As(eiz_z[j],eiz_w[j],length,n,y[j,k])
		A_2[j,k] = A_2s(eiz_z[j],eiz_w[j],length,n,y_2[j,k])# * np.cos(2.0*theta)  		
		A[0] = (1./2)*A[0]  		
		A_2[0] = (1./2)*A_2[0]  		
	sum_A = np.sum(A,axis=0)
	#print 'sum of A0 = ', k,j,sum_A
	sum_A_2 = np.sum(A_2,axis=0)
	#print 'sum of A2 = ', k,j,sum_A_2
	#print '----'
	#print 'shape sum_A_2 = ', np.shape(sum_A_2)
	#sys.exit()
for k,length in enumerate(ls):
	for i, theta in enumerate(thetas):
		EL[k] = 1./(length*length*length*length)
		G_l_t_dt[k,i] = (1.602e-19 / 4.11e-21) * (1./32) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A[k] + sum_A_2[k]* np.cos(2.0*theta) )/(2.0*np.sin(theta))# (1e21)*
np.savetxt('G_srw.txt',G_l_t_dt)

pl.figure()
pl.plot(ns,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
pl.plot(ns,eiz_y, color = 'g', label = r'$\varepsilon_{\hat{y}}(i\zeta_{N})$')
pl.plot(ns,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{N})$')
#pl.plot(ns,eiz_w, color = 'c', label = r'$\varepsilon_{vac}(i\zeta_{N})$')
pl.plot(ns,eiz_w, color = 'c', label = r'$\varepsilon_{water}(i\zeta_{N})$')
pl.xlabel(r'$N$', size = 20)
pl.ylabel(r'$\varepsilon(i\zeta)$', size = 20)
pl.legend(loc = 'best')
#pl.title(r'$\mathrm{CG-10\, DNA}$', size = 20)
pl.axis([0,500,0.9,2.6])
pl.savefig('plots/skew_ret_water/eiz.pdf' )
show()

pl.figure()
pl.plot(ns,aiz, color = 'b')#, label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
pl.xlabel(r'$N$', size = 20)
pl.ylabel(r'$a_{1,2}(i\zeta_{N})$', size = 20)
pl.legend(loc = 'best')
#pl.title(r'$\mathrm{Anisotropy \,Metric}$', size = 20)
#pl.axis([0,500,0.9,2.6])
pl.savefig('plots/skew_ret_water/aiz.pdf' )
show()

pl.figure()
pl.loglog(ls,(kb*T/32)*sum_A,'b-', label = r'$\mathcal{A^{(0)}}$')
pl.loglog(ls,(kb*T/32)*sum_A_2,'b--', label = r'$\mathcal{A^{(2)}}$')
pl.xlabel(r'$\mathrm{separation}\,\ell\,\,\,\rm{[m]}$', size = 20)
pl.ylabel(r'$\mathrm{\mathcal{A^{(0)},\,\,A^{(2)}}}$', size = 20)
#pl.title(r'$\mathrm{Hamaker \, coeff.s \,:\,skewed,\,retarded,\,water}$', size = 20)
pl.legend(loc = 'upper right')
pl.axis([1e-9,1e-6,1e-24,1e-19])
pl.savefig('plots/skew_ret_water/skew_ret_water_A0_A2.pdf')
show()

ls4 = 1e9*ls[ 2]#2]
ls5 = 1e9*ls[12]#4]
ls6 = 1e9*ls[22]#6]
ls1 = 1e9*ls[32]#8]
ls2 = 1e9*ls[42]#12]
ls3 = 1e9*ls[52]#16]

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
#pl.semilogy(thetas, G_l_t_dt)
ax.semilogy(thetas, G_l_t_dt[ 2,:], label = r'$\ell$ = %1.2f nm' %ls4)
ax.semilogy(thetas, G_l_t_dt[12,:], label = r'$\ell$ = %1.2f nm' %ls5)
ax.semilogy(thetas, G_l_t_dt[22,:], label = r'$\ell$ = %1.2f nm' %ls6)
ax.semilogy(thetas, G_l_t_dt[32,:], label = r'$\ell$ = %1.2f nm' %ls1)
ax.semilogy(thetas, G_l_t_dt[42,:], label = r'$\ell$ = %1.2f nm' %ls2)
ax.semilogy(thetas, G_l_t_dt[52,:], label = r'$\ell$ = %1.2f nm' %ls3)
#ax.semilogy(0,0,'', label = r'$G_\theta = cos(2\theta)/2sin(\theta)$')
pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
pl.axis([0,1.7,1e-10,1.0])
#pl.axis([0,1.7,1e-3,1e4])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,retarded,\,water}$', size = 20)
pl.legend(loc = 'lower left')
#pl.savefig('plots/skew_ret_water/skew_ret_water_G_vs_theta.pdf')
#show()
pl.savefig('plots/skew_ret_water/G_vs_theta_fixed_l.pdf')
show()

pl.figure()
#pl.loglog(ls, G_l_t_dt)#, label = labels[i])
pl.loglog(ls, G_l_t_dt[:,3], label = r'$\theta = \pi/4$')
pl.loglog(ls, G_l_t_dt[:,4], label = r'$\theta = \pi/3$')
pl.loglog(ls, G_l_t_dt[:,6], label = r'$\theta = \pi/2$')
pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
pl.axis([1.0e-9, 1.0e-6,1e-16,1e3])
#pl.axis([1.0e-9, 1.0e-6,1e-3,1e3])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,retarded,\,water}$', size = 20)
pl.legend(loc = 'best')
pl.savefig('plots/skew_ret_water/skew_ret_water_G_vs_l.pdf')
show()

#pl.figure()
##pl.semilogy(thetas, G_l_t_dt[12,:], label = r'$\ell$\,=\,4.0\,nm')
#pl.plot(thetas, G_l_t_dt[ 8,:], label = r'$\ell$ = %1.2f nm' %ls1)
#pl.plot(thetas, G_l_t_dt[12,:], label = r'$\ell$ = %1.2f nm' %ls2)
#pl.plot(thetas, G_l_t_dt[16,:], label = r'$\ell$ = %1.2f nm' %ls3)
#pl.plot(thetas, np.cos(2.*thetas)/(2.*np.sin(thetas)), label = r'$y = cos(2\theta)/2sin(\theta)$')
#pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
#pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
#pl.axis([0,6.28,-5,5])
##pl.axis([(1./25)*np.pi,(3./4)*np.pi,105,135])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.savefig('plots/skew_ret_water/G_vs_theta_fixed_l.pdf')
#show()
#ts1 = 1e9*1s[8]
#ts2 = 1e9*1s[12]
#ts3 = 1e9*1s[16]

#fig = pl.figure()
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
##pl.semilogy(thetas, G_l_t_dt)
#pl.semilogy(thetas, G_l_t_dt[ 8,:], label = r'$\ell$ = %1.2f nm' %ls1)
#pl.semilogy(thetas, G_l_t_dt[12,:], label = r'$\ell$ = %1.2f nm' %ls2)
#pl.semilogy(thetas, G_l_t_dt[16,:], label = r'$\ell$ = %1.2f nm' %ls3)
#pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
#pl.ylabel(r'$G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([(1./25)*np.pi,(3./4)*np.pi,105,135])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.savefig('plots/skew_ret_water/skew_ret_water_G_vs_theta.pdf')
#show()

#fig = pl.figure()
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#ax.plot(x_peak  , y_peak,   color = 'g',label = r'Synthesized')# $\epsilon$"($\omega$)', linestyle='-')
#ax.plot(x_nopeak, y_nopeak, color = 'b',label = r'Ab initio')# $\epsilon$"($\omega$)',   linestyle='-')
##ax.legend(loc = 'best')
##pl.title(r'$\epsilon$"($\omega$)  Ab Initio and Synthisized')
#pl.xlabel(r'$\omega$', size = 22)
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
#
#ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
#ax_inset.plot(x_nopeak,diff_eps,color='r',label=r'$\epsilon$"($\omega)_{synthesized}$-$\epsilon$"($\omega)_{ab\,initio}$')
#ax_inset.plot(x_nopeak,listofzeros,color = 'k', label=r'$\delta$$\epsilon$"($\omega$) = 0')
##pl.title(r'Difference $\epsilon$"($\omega$)', size = 'small')
#pl.tick_params(labelsize = 'small')
#pl.xlabel(r'$\omega$', size = 14)
#pl.ylabel(r'$\delta$$\epsilon$"($\omega$)', size = 14)
##pp.savefig()
##pl.savefig('130807_plots/fig2a_eps2_deps2.pdf')
##pl.show()


###ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
####pl.semilogy(thetas, G_l_t_dt[12,:], label = r'$\ell$\,=\,4.0\,nm')
###ax_inset.plot(thetas, G_l_t_dt[ 8,:], label = r'$\ell$ = %1.2f nm' %ls1)
###ax_inset.plot(thetas, G_l_t_dt[12,:], label = r'$\ell$ = %1.2f nm' %ls2)
###ax_inset.plot(thetas, G_l_t_dt[16,:], label = r'$\ell$ = %1.2f nm' %ls3)
###ax_inset.plot(thetas, np.cos(2.*thetas)/(2.*np.sin(thetas)))#, label = r'$G_\theta = cos(2\theta)/2sin(\theta)$')
###pl.tick_params(labelsize = 'small')
###pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 10)
###pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 10)
###pl.axis([0,2.89,-5,5])
#pl.axis([(1./25)*np.pi,(3./4)*np.pi,105,135])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'best')

#pl.figure()
#pl.loglog(thetas, G_l_t_dt)#, label = labels_l[k])
#pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
#pl.ylabel(r'$G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([(1./25)*np.pi,(3./4)*np.pi,105,135])
#pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,retarded,\,water}$', size = 20)
##pl.legend(loc = 'best')
#pl.savefig('plots/skew_ret_water/skew_ret_water_G_vs_theta.pdf')
#show()
