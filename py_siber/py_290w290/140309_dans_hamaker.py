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
#pp = PdfPages('plots/skew_ret_water/skew_ret_water.pdf')

#eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') #perpendicular, radial
#eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # parallel,axial
#eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # water as intervening medium
##eiz_w[0] = eiz_w[1] #NOTE: there is a jump from first val down to second val

gh_len0, gh_A0 = np.loadtxt('data/29-w-29-500nmA0.txt',unpack=True, usecols = [0,1]) # water as intervening medium
gh_len2, gh_A2 = np.loadtxt('data/29-w-29-500nmA2.txt',unpack=True, usecols = [0,1]) # water as intervening medium
                                                       #
x65x00_A0, y65y00_A0 = np.loadtxt('data/65-W-65-PARA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x65x00_A2, y65y00_A2 = np.loadtxt('data/65-W-65-PARA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                        
x91x00_A0, y91y00_A0 = np.loadtxt('data/91-W-91-PARA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x91x00_A2, y91y00_A2 = np.loadtxt('data/91-W-91-PARA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                        
x93x00_A0, y93y00_A0 = np.loadtxt('data/93-W-93-PARA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x93x00_A2, y93y00_A2 = np.loadtxt('data/93-W-93-PARA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                      
X290x00_A0, y290y00_A0 = np.loadtxt('data/290-W-290-PARA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
X290x00_A2, y290y00_A2 = np.loadtxt('data/290-W-290-PARA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                      
x65x90_A0,   y65y90_A0 =   np.loadtxt('data/65-W-65-PERPA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x65x90_A2,   y65y90_A2 =   np.loadtxt('data/65-W-65-PERPA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                        
x91x90_A0,   y91y90_A0 =   np.loadtxt('data/91-W-91-PERPA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x91x90_A2,   y91y90_A2 =   np.loadtxt('data/91-W-91-PERPA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                        
x93x90_A0,   y93y90_A0 =   np.loadtxt('data/93-W-93-PERPA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x93x90_A2,   y93y90_A2 =   np.loadtxt('data/93-W-93-PERPA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                        
X290x90_A0, y290y90_A0 = np.loadtxt('data/290-W-290-PERPA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
X290x90_A2, y290y90_A2 = np.loadtxt('data/290-W-290-PERPA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
                                                      
x90_A0, y90_A0 = np.loadtxt('data/290-W-290-PERPA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x90_A2, y90_A2 = np.loadtxt('data/290-W-290-PERPA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium

pl.figure()
#pl.loglog(gh_len0, gh_A0,'k:', marker = 'o', markevery = 3, label = r'GH $\mathcal{A^{(0)}}(\ell)$')
#pl.loglog(gh_len2, gh_A2,'k:', marker = 'o', markevery = 3,label = r'GH $\mathcal{A^{(2)}}(\ell)$')

pl.loglog( x65x00_A0,  y65y00_A0,'b:' , linewidth = 3, marker = '+', markevery = 3, label = r'$\mathcal{A^{(0)}}_{65w65}(\ell, \theta = 0)$')
pl.loglog( x65x00_A2,  y65y00_A2,'b--', linewidth = 1, marker = '+', markevery = 3, label = r'$\mathcal{A^{(2)}}_{65w65}(\ell, \theta = 0)$')
                                                                            
pl.loglog( x91x00_A0,  y91y00_A0,'g:' , linewidth = 3, marker = 'x', markevery = 3, label = r'$\mathcal{A^{(0)}}_{91w91}(\ell, \theta = 0)$')
pl.loglog( x91x00_A2,  y91y00_A2,'g--', linewidth = 1, marker = 'x', markevery = 3, label = r'$\mathcal{A^{(2)}}_{91w91}(\ell, \theta = 0)$')
                                                                            
pl.loglog( x93x00_A0,  y93y00_A0,'r:' , linewidth = 3, marker = '.', markevery = 3, label = r'$\mathcal{A^{(0)}}_{93w93}(\ell, \theta = 0)$')
pl.loglog( x93x00_A2,  y93y00_A2,'r--', linewidth = 1, marker = '.', markevery = 3, label = r'$\mathcal{A^{(2)}}_{93w93}(\ell, \theta = 0)$')
                       
pl.loglog(X290x00_A0, y290y00_A0,'c:' , linewidth = 3, marker = '|', markevery = 3, label = r'$\mathcal{A^{(0)}}_{290w290}(\ell, \theta = 0)$')
pl.loglog(X290x00_A2, y290y00_A2,'c--', linewidth = 1, marker = '|', markevery = 3, label = r'$\mathcal{A^{(2)}}_{290w290}(\ell, \theta = 0)$')


pl.loglog( x65x90_A0,  y65y90_A0,'k:', marker = '+', markevery = 10, label = r'$\mathcal{A^{(0)}}_{65w65}(\ell, \theta =\pi/2)$') 
pl.loglog( x65x90_A2,  y65y90_A2,'k:', marker = '+', markevery = 10, label = r'$\mathcal{A^{(2)}}_{65w65}(\ell, \theta =\pi/2)$')  
                                                                                                                               
pl.loglog( x91x90_A0,  y91y90_A0,'k:', marker = 'x', markevery = 10, label = r'$\mathcal{A^{(0)}}_{91w91}(\ell, \theta =\pi/2)$')  
pl.loglog( x91x90_A2,  y91y90_A2,'k:', marker = 'x', markevery = 10, label = r'$\mathcal{A^{(2)}}_{91w91}(\ell, \theta =\pi/2)$')  
                                                                                                                               
pl.loglog( x93x90_A0,  y93y90_A0,'k:', marker = '.', markevery = 10, label = r'$\mathcal{A^{(0)}}_{93w93}(\ell, \theta =\pi/2)$')  
pl.loglog( x93x90_A2,  y93y90_A2,'k:', marker = '.', markevery = 10, label = r'$\mathcal{A^{(2)}}_{93w93}(\ell, \theta =\pi/2)$')  
                                                                                                                               
pl.loglog(X290x90_A0, y290y90_A0,'k:', marker = '|', markevery = 10, label = r'$\mathcal{A^{(0)}}_{290w290}(\ell, \theta =\pi/2)$')
pl.loglog(X290x90_A2, y290y90_A2,'k:', marker = '|', markevery = 10, label = r'$\mathcal{A^{(2)}}_{290w290}(\ell, \theta =\pi/2)$')
pl.axis([1e0,0.2e4,1e-1,1e2])

pl.xlabel(r'$\ell\,\,\,\rm{[nm]}$', size = 20)
pl.ylabel(r'$\mathrm{\mathcal{A^{(0)},\,\,A^{(2)}}}\,\,\,\rm{[zJ]}$', size = 20)
pl.legend(loc = 'center right')
pl.title('Gecko Hamaker coefficients')
pl.show()

