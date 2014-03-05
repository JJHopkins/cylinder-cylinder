#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from scipy import integrate
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('130805_Hopkins_ab_synth_eps2_eiz_normed_y_and_df_DF.pdf')

pl.close('all')
#x_t ,y_t = np.loadtxt('data/CG10e2t.csv',delimiter=',',unpack=True, usecols = [0,1])
#x_x ,y_x= np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])
#x_y ,y_y = np.loadtxt('data/CG10e2y.csv',delimiter=',',unpack=True, usecols = [0,1])
#x_z ,y_z = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])
#
#eiz_x = np.loadtxt('data/eiz_x_output.txt',unpack=True, usecols = [0])
#eiz_z = np.loadtxt('data/eiz_z_output.txt',unpack=True, usecols = [0])

x_nopeak,y_nopeak = np.loadtxt('data/CG10e2z.csv',delimiter=',',unpack=True, usecols = [0,1])
x_peak,y_peak = np.loadtxt('data/CG10e2x.csv',delimiter=',',unpack=True, usecols = [0,1])

eiz_nopeak = np.loadtxt('data/eiz_z_output.txt',unpack=True, usecols = [0])
eiz_peak = np.loadtxt('data/eiz_x_output.txt',unpack=True, usecols = [0])

coeff = 0.159 # in eV  #(2.41*1e14) # in rad/s
n = arange(0,250)
z = n * coeff
z_y = linspace(0,250,250)
A = 4.0 #[1,4,6]
B = 7.0 #[1,4,6]
C = 10.0 #[1,4,6]

sum_Li2 = numpy.empty(len(z))
dF = numpy.empty(len(z))

sum_Li2_A_dep = numpy.empty(len(z_y))

yA_dep = numpy.empty(len(z_y))
norm_yA_dep = numpy.empty(len(z_y))

A_dep = linspace(1,51,100)
sum_dF =  0.0
sum_F =  0.0
sum_F_p =  0.0
diff_sum_F =  0.0
max_yA = numpy.empty(len(A_dep))
sum_Li3 = numpy.empty(len(z))
sum_Li3_p = numpy.empty(len(z))
F_3 = numpy.empty(len(z))
F_3_p = numpy.empty(len(z))

for j in range(len(z)):
    sum_Li2[j] = 0.0
    for m in arange(101) + 1:
        sum_Li2[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**m)/(m**2)
    dF[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    sum_Li3[j] = 0.0
    sum_Li3_p[j] = 0.0
    
    for h in arange(101)+ 1:
        sum_Li3[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**h)/(h**3)
        sum_Li3_p[j] += ((((eiz_peak[j] -1)/(eiz_peak[j] +1))**2)**h)/(h**3)
sum_dF = sum(dF)
sum_F = sum(sum_Li3)
sum_F_p = sum(sum_Li3_p)
diff_sum_F = sum_F- sum_F_p
perc_diff_df = 0.0
perc_diff_df = (sum_dF-diff_sum_F)/diff_sum_F

print ('Total dF   = %s ' % sum_dF)
print ('Total F_np = %s ' % sum_F)
print ('Total F_p  = %s ' % sum_F_p)
print ('Total DF   = %s ' % diff_sum_F)
print ('Agreement: (dF-DF)/DF = %s ' % perc_diff_df)

dF[0] = (1./2)*dF[0]   
sum_Li3[0] = (1./2)*sum_Li3[0]
sum_Li3_p[0] = (1./2)*sum_Li3_p[0]

F_3 = -sum_Li3
F_3_p = -sum_Li3_p
diff_F = (-F_3+F_3_p)
divide = dF/diff_F 

## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
diff_eps = y_peak - y_nopeak
diff_eiz = eiz_peak - eiz_nopeak
listofzeros = numpy.zeros(len(x_nopeak)) # plot line for y = 0

## PLOTS
#-------------------------------------------------------------
#pl.figure()
#pl.plot(z_y,yA/yA[0], color = 'b', label = 'A = 1.0')
#pl.plot(z_y,yB/yB[0], color = 'g', label = 'B = 4.0')
#pl.plot(z_y,yC/yC[0], color = 'r', label = 'C = 6.0')
#pl.title(r'Function $y(\xi)$ from $eqn\,6$')
#pl.tick_params(labelsize = 'small')
#pl.xlabel(r'$0<\frac{\xi}{\omega_0}<2$', size = 'x-large')
#pl.ylabel(r'$y(\xi)$', size = 'x-large')# \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$')
#pl.legend(loc = 'best')
#pl.axis([0.0, 2.2,0.0,1.2])
#pp.savefig()
#pl.show()
##############################################################
fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(x_peak  , y_peak,   color = 'g',label = r'Synthesized')# $\epsilon$"($\omega$)', linestyle='-')
ax.plot(x_nopeak, y_nopeak, color = 'b',label = r'Ab initio')# $\epsilon$"($\omega$)',   linestyle='-')
#ax.legend(loc = 'best')
#pl.title(r'$\epsilon$"($\omega$)  Ab Initio and Synthisized')
pl.xlabel(r'$\omega$', size = 22)
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)

ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
ax_inset.plot(x_nopeak,diff_eps,color='r',label=r'$\epsilon$"($\omega)_{synthesized}$-$\epsilon$"($\omega)_{ab\,initio}$')
ax_inset.plot(x_nopeak,listofzeros,color = 'k', label=r'$\delta$$\epsilon$"($\omega$) = 0')
#pl.title(r'Difference $\epsilon$"($\omega$)', size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$\omega$', size = 14)
pl.ylabel(r'$\delta$$\epsilon$"($\omega$)', size = 14)
#pp.savefig()
#pl.savefig('130807_plots/fig2a_eps2_deps2.pdf')
#pl.show()
################################################################

# ab synth and diff in eV
fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(  x_peak,  y_peak,   color = 'g',label = r'Synthesized')# $\epsilon$"($\omega$)', linestyle='-')
ax.plot(x_nopeak,y_nopeak, color = 'b',label = r'Ab initio')# $\epsilon$"($\omega$)',   linestyle='-')
#ax.legend(loc = 'best')
#pl.title(r'$\epsilon$"($\omega$)  Ab Initio and Synthisized')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)

ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
ax_inset.plot(x_nopeak,diff_eps,color='r',label=r'$\epsilon$"($\omega)_{synthesized}$-$\epsilon$"($\omega)_{ab\,initio}$')
ax_inset.plot(x_nopeak,listofzeros,color = 'k', label=r'$\delta$$\epsilon$"($\omega$) = 0')
#pl.title(r'Difference $\epsilon$"($\omega$)', size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 14)
pl.ylabel(r'$\delta$$\epsilon$"($\omega$)', size = 14)
#pp.savefig()
#pl.savefig('plots/130814_Hopkins_eV_ab_initio_synth_eps2_delta_eps2.pdf')
####pl.savefig('Fig2a.pdf')
####pl.savefig('plots/130911_Hopkins_ab_initio_synth_eps2_delta_eps2.png', dpi = 300)
#pl.show()
#################################################################

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(n, eiz_peak  , color = 'g', label = r'Synthesized')#$ \epsilon (i \zeta ) _{peak}$')
ax.plot(n, eiz_nopeak, color = 'b', label = r'Ab initio')#$ \epsilon (i \zeta )_{no \,peak} $')
#ax.legend(loc = 'best')
#pl.title('$\epsilon$(i$\zeta$) for synthesized and ab initio data')
pl.xlabel(r'$N$', size = 22)
pl.ylabel(r'$ \epsilon (i \zeta_{N} ) $', size = 21)
#pl.axis([0.0, 150,1.3,2.4])

ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
ax_inset.plot(n,diff_eiz,color= 'r',linestyle='-')
#pl.title('Difference $\epsilon$(i$\zeta$)',size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$N$',size = 14)#(r'$\zeta$')
pl.ylabel(r'$ \delta \epsilon (i  \zeta_{N} )$',size = 14)
#pl.axis([0.0, 150,0.0,0.1])
#pp.savefig()
####pl.savefig('Fig2b.pdf')
#####pl.savefig('plots/130911_Hopkins_ab_initio_synth_eiz_delta_eiz.png', dpi = 300)
#pl.savefig('130807_plots/fig2b_eiz_deiz.pdf')
#####################################################################
#pl.figure()
fig = pl.figure()#figsize=(6,4.5))
ax = fig.add_axes([0.14,0.1,0.8,0.8])
ax.plot(n,-F_3  , color = 'b', label = r'$\left( F/S \right)_{ab\,initio}$')
ax.plot(n,-F_3_p, color = 'g', linestyle='-',label=r'$\left( F/S \right) _{synthesized}$')
#pl.title(r'vdW Free Energy as a function of n$^{th}$ Matsubara fequency')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$N$', size = 21)
pl.ylabel(r'$\mathcal{A}\mathrm{(i\zeta_{N})/k _{B} T}$', size = 21)
#pl.ylabel(r'$\left( \mathrm{F/S} \right) \,\, \frac{ 12\pi D^{2} } {k _{B} T}$', size = 21)
#pl.legend(loc='best')
#pp.savefig()
#pl.show()
#pl.axis([0.02, 150,0, 0.22])

#pl.figure()
ax_inset = fig.add_axes([0.57,0.51,0.36,0.36])
ax_inset.plot(n,-dF, color = 'r', label = r'$\delta F$')
#ax_inset.plot(n,diff_F, color = 'k', linestyle = '--', label = r'$\Delta F$')
#pl.title(r'vdW Free Energy change as a function of n$^{th}$ Matsubara fequency')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$N$', size = 14)
pl.ylabel(r'$\delta\mathcal{A}\mathrm{(i\zeta_{N})/k _{B} T}$', size = 14)
#pl.ylabel(r'$\delta \mathrm{ F}\,\, \frac{ 12\pi D^{2} }{k_{B} T}$')
#pl.legend(loc = 'best')
#pl.axis([0.0, 150,0.0,0.016])
#pp.savefig()
#pl.savefig('130807_plots/fig3_F_dF.pdf')
####pl.savefig('Fig3.pdf')
###pl.savefig('plots/130911_Hopkins_ab_initio_synth_A_deltaA.png', dpi = 300)
pl.show()

#pp.close()
pl.close()






