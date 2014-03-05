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
pp = PdfPages('plots/par_NR_water/par_NR_water.pdf')

eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') #perpendicular, radial
eiz_y = np.loadtxt('data/eiz_y_output_eV.txt')
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # parallel,axial
#eiz_w = 1.0 + np.zeros(len(eiz_z))
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # water as intervening medium

eiz_w[0] = eiz_w[1] #NOTE: there is a jump from first val down to second val

r_1 = 0.5e-9
r_2 = 0.5e-9
c = 2.99e8 # in m/s
# at RT, 1 kT = 4.11e-21 J
# 1 eV = 1.602e-19 J = 0.016 zJ
# h_bar_eV = 6.5821e-16 eVs
T = 297 
# h_bar = 1. #1.0546e-34 #in Js
#kb = 8.6173e-5  # in eV/K
kb  = 1.3807e-23 # in J/K
# NOTES:
# z_n_eV = (2*pi*kT/h_bar)n 
#	= (0.159 eV) / (6.5821e-16 eVs) 
#	= n*2.411e14 rad/s
# z_n_J = (2*pi*kT/h_bar)n 
#	= (1.3807e-23 J/K) / (1.0546e-34 Js))*n
#	= n*2.411e14 rad/s
#coeff = 0.159 # in eV w/o 1/h_bar
coeff = 2.411e14 # in rad/s

ns = np.arange(0.,500.)
z = ns * coeff
ls = np.linspace(1.0e-9, 7.0e-9, 25)

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
def ys(a):
	term1 =	3.0 * 5.*(a + a)
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

A   = np.zeros(shape=(len(ns),len(ls)))
A_2 = np.zeros(shape=(len(ns),len(ls)))
aiz = []
sum_A = np.empty(len(ls))
sum_A_2 = np.empty(len(ls))
EL = np.zeros(len(ls))
G_l_t_dt = np.empty(len(ls))

aiz = Aiz(eiz_x,eiz_z, eiz_w) # of length = len(ns)

#sum_A = np.empty(len(ls))
#sum_A_2 = np.empty(len(ls))
#for j,n in enumerate(ns):
#	# Integral:
#	y[j]   = ys(aiz[j])
#	y_2[j] = y_2s(aiz[j])
#	#print 'dt Integral   y = ',i,k,j, y
#	#print 'dt Integral y_2 = ',i,k,j, y_2
#	#print '----'
#	#print 'N terms for A0 = '  , As(eiz_z[j],eiz_w[j],length,n,y)
#	#print 'N terms for A2 = ', A_2s(eiz_z[j],eiz_w[j],length,n,y_2)
#	#print '----'
#	A[j]   = As(eiz_z[j],eiz_w[j],y[j])
#	A_2[j] = A_2s(eiz_z[j],eiz_w[j],y_2[j])# * np.cos(2.0*theta)  		
#sum_A = np.sum(A,axis=0)
#print 'sum of A0 = ', j,sum_A
#sum_A_2 = np.sum(A_2,axis=0)
#print 'sum of A2 = ', j,sum_A_2
#print '----'
#for k,length in enumerate(ls):
#	for i, theta in enumerate(thetas):
#		EL[k] = 1./(length*length*length*length)
#		#print '-------------------------------------i,theta = ', (i,theta)
#		G_l_t_dt[k,i] = (1.602e-19 / 4.11e-21) * (1./32) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A + sum_A_2* np.cos(2.0*theta) )/(2.0*np.sin(theta))# (1e21)*
#		#G_l_t_dt[k,i] = (1.602e-19 / 4.11e-21) * (1./32) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A[k] + sum_A_2[k]* np.cos(2.0*theta) )/(2.0*np.sin(theta))# (1e21)*
#		#labels = 'theta = %.3f' %(theta) 
#		#labels_l = 'distance = %1.2f nm' %(length*1e9) 
#        	#print 'theta = %.1f, length = %.1f, G = %s' %(i,k,G_l_t_dt[k,i]) 
### 1 eV = 1.602e-19 J = 0.016 zJ
### h_bar_eV = 6.5821e-16 eVs
for k,length in enumerate(ls):
	for j,n in enumerate(ns):
		#print "on n=%d of %d"%(j,len(ns))
		# Integrand:
		y   = ys(aiz[j])
		y_2 = y_2s(aiz[j])
		A[j,k]   = As(eiz_z[j],eiz_w[j],y)
		A_2[j,k] = A_2s(eiz_z[j],eiz_w[j],y_2)  		
		A[0,k] = (1./2)*A[0,k]  		
		A_2[0,k] = (1./2)*A_2[0,k]  		
	sum_A = np.sum(A,axis=0)
	#print 'shape sum_A = ', np.shape(sum_A)
	sum_A_2 = np.sum(A_2,axis=0)
	#print 'shape sum_A_2 = ', np.shape(sum_A_2)
	#sys.exit()
	EL[k] = 1./(length*length*length*length*length)
	G_l_t_dt[k] = kb * T * (9./(32.*64.)) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A[k] + sum_A_2[k])
	#print 'theta = %.1f, length = %.1f, G = %s' %(i,k,G_l_t_dt[k,i]) 

pl.figure()
pl.plot(ns,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
pl.plot(ns,eiz_y, color = 'g', label = r'$\varepsilon_{\hat{y}}(i\zeta_{N})$')
pl.plot(ns,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{N})$')
#pl.plot(ns,eiz_w, color = 'c', label = r'$\varepsilon_{vac}(i\zeta_{N})$')
pl.plot(ns,eiz_w, color = 'c', label = r'$\varepsilon_{water}(i\zeta_{N})$')
pl.xlabel(r'$N$', size = 20)
pl.ylabel(r'$\varepsilon(i\zeta)$', size = 20)
#pl.title('')
pl.legend(loc = 'best')
pl.title(r'$\mathrm{CG-10\, DNA}$', size = 20)
pl.axis([0,500,0.9,2.6])
pl.savefig('plots/par_NR_water/eiz.pdf' )
pl.show()
show()
#pl.figure()
#pl.loglog(ls,(kb*T/32)*sum_A,'b-', label = r'$\mathcal{A^{(0)}}$')
#pl.loglog(ls,(kb*T/32)*sum_A_2,'b--', label = r'$\mathcal{A^{(2)}}$')
#pl.xlabel(r'$\mathrm{separation}\,\ell\,\,\,\rm{[m]}$', size = 20)
#pl.ylabel(r'$\mathrm{\mathcal{-A^{(0)},\,\,-A^{(2)}}}$', size = 20)
#pl.title(r'$\mathrm{Hamaker \, coeff.s \,:\,skewed,\,non-ret,\,water}$', size = 20)
##pl.title(r'$\mathrm{Hamaker \, coeff.s \,\mathcal{-A^{(0)},\,\,-A^{(2)}}\,:\,skewed,\,non-ret,\,water}$', size = 20)
##pl.title(r'$\mathrm{Hamaker \, coefficients \,\mathcal{A_{0},A_{2}}}\, skewed,\,non-ret,\,water$', size = 20)
#pl.legend(loc = 'lower right')
#pl.axis([1e-9,1e-8,1e-24,1e-18])
#pl.savefig('plots/par_NR_water/A0_A2.pdf')
#show()


pl.figure()
pl.plot(ls,sum_A+sum_A_2,'b-')
pl.xlabel(r'$\mathrm{separation}\,\it{l}\,\,\,\rm{[nm]}$', size = 20)
pl.ylabel(r'$\mathrm{\mathcal{A_{0},A_{2}}\,[J]}$', size = 20)
pl.title(r'$\mathrm{Hamaker \, coefficients \,\mathcal{A_{0},A_{2}}}\, parallel,\,retarded,\,water$', size = 20)
show()

pl.figure()
pl.loglog(ls, G_l_t_dt)
pl.xlabel(r'$Separation\,\,[m]$', size = 24)
pl.ylabel(r'$g(l)\,\,[J]$', size = 24)
#pl.axis([1.5e-9, 6.5e-8,100,145])
pl.title('g as a function of separation')
show()

pl.figure()
pl.plot(ls, G_l_t_dt)
pl.xlabel(r'$Separation\,\,[m]$', size = 24)
pl.ylabel(r'$g(l)\,\,[J]$', size = 24)
#pl.axis([1.5e-9, 6.5e-8,100,145])
pl.title('g as a function of separation')
show()
pl.figure()
#pl.loglog(thetas, G_l_t_dt)#, label = labels_l[k])
pl.semilogy(thetas, G_l_t_dt)
pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
pl.ylabel(r'$G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
#pl.axis([(1./25)*np.pi,(3./4)*np.pi,105,135])
#pl.title('G as a function of angle')
pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,non-ret,\,water}$', size = 20)
#pl.legend(loc = 'best')
pl.savefig('plots/skew_NR_water/semilog_G_vs_theta.pdf')
show()

pl.figure()
pl.loglog(ls, G_l_t_dt)#, label = labels[i])
pl.xlabel(r'$Separation,\,\ell\,\,\mathrm{[m]}$', size = 20)
pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
#pl.axis([1.5e-9, 6.5e-8,100,145])
#pl.title('G as a function of separation')
pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,non-ret,\,water}$', size = 20)
#pl.legend(loc = 'best')
pl.savefig('plots/skew_NR_water/G_vs_l.pdf')
show()

#G_l_t_dt[G_l_t_dt>300]= np.nan #NOTE: remove me later
#G_l_t_dt[G_l_t_dt<200e-25]= np.nan #NOTE: remove me later
# CONTOUR PLOT:
X,Y = np.meshgrid(thetas, ls)
pl.figure()
pl.contourf(X, Y, np.log(G_l_t_dt), 25)#, cmap = cm.hot())
CS = pl.contour(X,Y,np.log(G_l_t_dt))#, levels = np.linspace(1e-1,1e10,10))
pl.clabel(CS, inline =1,fmt = '%1.5f', fontsize = 18,color = 'k')#, manual = man_loc)
pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
pl.ylabel(r'$Separation\,\,\mathrm{[m]}$', size = 20)
pl.title(r'$\mathrm{-Log(G),\,non-ret,\,skewed\,cyls\,in\,water}$', size = 20)#uas a function of separation and angle')
cbar = pl.colorbar(CS, shrink = 0.8, extend = 'both')
cbar.ax.set_ylabel(r'$-Log(G(\mathcal{\ell},\theta))\,\,[k_{B}T]$', size = 14)
cbar.add_lines(CS)
##pl.axis([0,1.0,0,1.0])
#pl.grid()
pl.savefig('plots/skew_NR_water/logG_contour.pdf')
show()

fig = pl.figure()
ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
surf = ax.plot_surface(X,Y, G_l_t_dt, rstride = 1, cstride = 1,alpha = 0.2, linewidth = 0.3)#edgecolor = 'none',antialiased = True, shade = False, norm = norm, linewidth = 0.3)
#surf = ax.plot_surface(X,Y, G_l_t_dt, rstride = 20, cstride = 20,alpha = 0.2)#, cmap = cm.gnuplot, linewidth = 0.5)#gray)#coolwarm)#bone)#hot, linewidth = 0.01, antialiased = True, shade = False)# True)#, cmap = hot()
#colorbar(surf)
#cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
#cset = ax.contour(X,Y,h, zdir = 'z', offset = 0, cmap = cm.jet)
#cset = ax.contour(X,Y,h, zdir = 'x', offset = 5, cmap = cm.jet)
#cset = ax.contourf(X,Y,h, zdir = 'y', offset = 6, cmap = cm.jet)# puts plot of max xi vs discrete r values at r=0 plane
#ax.view_init(elev = 19, azim = -112)
#zlabel(r'$\xi/\omega_{0}$', size = 21)
#ylabel(r'$r$', size = 24)
#xlabel(r'$(\epsilon(0) -1)$', size = 24)
#text = Axes.text(self, x, y, s, **kwargs)
#art3d.text_2d_to_3d(text, z, zdir)
#return text
#pl.text(6,0, 0, r'$\xi/\omega_{0}$',size = 21 ,rotation = 'horizontal')
#ax.text(r'$\xi/\omega_{0}$',6,0, 0, size = 21 ,rotation = 'horizontal')
#ax.set_zlabel(r'$\xi/\omega_{0}$',size = 21 ,rotation = 'horizontal' )
ax.set_xlabel(r'$\epsilon(0)-1$', size = 21)
ax.set_ylabel(r'$r$', size = 22)
show()
#pp.savefig()




