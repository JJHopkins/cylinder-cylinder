#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130805_raz120_Hopkins_contour_xiMax_of_r_A.pdf')


z_y = linspace(  0.0000000001,1.2,500)# z_y is z_n/w_0 ~= e16/e16
A_dep = linspace(0.0000000001, 12,500)#(0.000001,10,20)
rho = linspace(  0.0000000001,  6,500)

sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))

#A_dep, rho = meshgrid(A_dep, rho)
height=zeros(shape=(len(A_dep),len(rho)))
for r in range(len(A_dep)):
    print >> sys.stdout, "on r=%d of %d"%(r,len(A_dep))
    for s in range(len(rho)):    

        for t in range(len(z_y)):
            sum_Li2_A_dep[t] = 0.0
	    #A_dep, rho = meshgrid(A_dep, rho)
	 
            ###for u in arange(1,50):	
	    U=arange(1,101,1)
		###sum_Li2_A_dep[t] += ((((A_dep[r]/(1.+(z_y[t]*rho[s])**2))/(2. + (A_dep[r]/(1. + (z_y[t]*rho[s])**2))))**2)**u)/u**2 
	    sum_Li2_A_dep[t] = sum ( ((((A_dep[r]/(1.+(z_y[t]*rho[s])**2))/(2. + (A_dep[r]/(1. + (z_y[t]*rho[s])**2))))**2)**U)/U**2 )
                #sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
            #print sum_Li2_A_dep[t]
            #yA_dep[t] = sum_Li2_A_dep[t]*\
            #(A_dep[r]/(A_dep[r]+2.*z_y[p]**2+2))**(-1)*\
            #(1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2)   
            
            yA_dep[t] = (sum_Li2_A_dep[t])*\
            ((1. - ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2))))**2)/\
            ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2)))))*\
            (1./(1. + z_y[t]**2))/(1. + A_dep[r]/(1+(z_y[t]*rho[s])**2)) #  NB this line is different, RP wrote numerator is 1/(1+x^2)...?
          
        height[r,s] = z_y[yA_dep==max(yA_dep)] 

X,Y=meshgrid(A_dep,rho)
figure()
levels = [0.20,0.40,0.60,0.80,1.00]
contourf(X,Y,height,1000); cbar = colorbar(); #clim(0,1);
CS=contour(X,Y,height, levels, colors = 'k')
clabel(CS,fmt = '%2.1f',colors = 'k',fontsize = 18)
xlabel(r'$A$', size = 22)#\,=\,\frac{2}{\pi}C$', size = 20)
ylabel(r'$r$', size = 22)
cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
cbar.add_lines(CS)
#show()

#fig = figure()
#ax = fig.add_subplot(111, projection = '3d')
#ax = Axes3D(fig)

#Axes3D.plot_wireframe(X,Y,height, rstride = .01, cstride = .01)
#ax.contour(A_dep,rho,max_xA)

##import sys
##sys.exit()

        
##max_rhoA = numpy.empty(len(rho))
#max_xA = numpy.empty(len(A_dep))
##max_rhoA[s] = rho[yA_dep.argmax()] 
#max_xA[r]= z_y[yA_dep.argmax()]
##print max_rhoA
#print max_xA
#        fig = figure()
#        ax = fig.add_subplot(111, projection = '3d')
#ax = Axes3D(fig)

#        Axes3D.plot_wireframe(A_dep,rho,max_yA, rstride = .01, cstride = .01)
#ax.contour(A_dep,rho,max_xA)
#data = zip(A_dep, max_xA)
#savetxt('output/130729_Hopkins_ymax_output.txt', data)

#pl.figure()
#pl.plot(A_dep, max_xA)  
##pl.plot(A_dep,max_yA/yA_dep[0], color = colors[r])
#pl.xlabel(r'$A\,=\,\frac{2}{\pi}C$', size = 'x-large')
#pl.ylabel(r'MAX$[y(\frac{\xi}{\omega_{0}})] $', size = 'x-large')
##pl.axis([0, 10.2, 0.0, 2.0])
pp.savefig()
#close()
pp.close()


