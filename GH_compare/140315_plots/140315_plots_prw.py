#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

A_py_par = r'$\mathcal{A}_{\parallel}\sf{[python]}$'
A_GH_par = r'$\mathcal{A}_{\parallel}\sf{[ G.H. ]}$'

x_ax = r'$\,\ell\,\,\,\rm{[nm]}$'
y_ax_par = r'$\mathrm{\mathcal{A}_\parallel(\ell)}\,\,\,\rm{[zJ]}$'

def title(cnt1,cnt2,orientation):
	return r'$\mathrm{[%s,%s]\,\,Hamaker\,coeff:\, %s \,in\,water,\,retarded}$'%(cnt1,cnt2,orientation)

def svfig(cnt1,cnt2,orientation):
	return 'plots/140316_%sw%s_HCs_%s_ret.pdf'%(cnt1,cnt2,orientation)


def plsvfig(x_p, y_p, A_p_lab, A_g_lab, x_lab,y_lab, Title_fnc,save_fnc):
	pl.figure()
	pl.loglog((1e9)*x_p, (1e21)*y_p,'b-', label = A_p_lab)
	pl.loglog((1e9)*x_p, (1e21)*y_p,'b-', label = A_g_lab)
	pl.xlabel(x_lab)
	pl.ylabel(y_lab)
	pl.title(Title_fnc)
	pl.legend(loc = 'best')
	pl.minorticks_on()
	pl.ticklabel_format(axis = 'both')
	pl.grid(which = 'both')
	pl.tick_params(which = 'both',labelright = True)
	pl.savefig(save_fnc)
	pl.show()

	return 0






#--- parallel ----------------------------------------------------------------
x_065_par = np.loadtxt('data/140316_65w65_prw_Ls.txt'   ) 
y_065_par = np.loadtxt('data/140316_65w65_prw_sum_A.txt') 
#
x_090_par = np.loadtxt('data/140316_90w90_prw_Ls.txt'   ) 
y_090_par = np.loadtxt('data/140316_90w90_prw_sum_A.txt') 
#
x_091_par = np.loadtxt('data/140316_91w91_prw_Ls.txt'   ) 
y_091_par = np.loadtxt('data/140316_91w91_prw_sum_A.txt') 
#
x_093_par = np.loadtxt('data/140316_93w93_prw_Ls.txt'   ) 
y_093_par = np.loadtxt('data/140316_93w93_prw_sum_A.txt') 
#
x_290_par = np.loadtxt('data/140316_290w290_prw_Ls.txt'   ) 
y_290_par = np.loadtxt('data/140316_290w290_prw_sum_A.txt') 

#--- perpendicular -----------------------------------------------------------
#x_065_per = np.loadtxt('data/140316_65w65_prw_Ls.txt'   ) 
#y_065_per = np.loadtxt('data/140316_65w65_prw_sum_A.txt') 
#
#x_090_per = np.loadtxt('data/140316_90w90_prw_Ls.txt'   ) 
#y_090_per = np.loadtxt('data/140316_90w90_prw_sum_A.txt') 
#
#x_091_per = np.loadtxt('data/140316_91w91_prw_Ls.txt'   ) 
#y_091_per = np.loadtxt('data/140316_91w91_prw_sum_A.txt') 
#
#x_093_per = np.loadtxt('data/140316_93w93_prw_Ls.txt'   ) 
#y_093_per = np.loadtxt('data/140316_93w93_prw_sum_A.txt') 
#
#x_290_per = np.loadtxt('data/140316_290w290_prw_Ls.txt'   ) 
#y_290_per = np.loadtxt('data/140316_290w290_prw_sum_A.txt') 


#--- 65w65 ------------------------------------------------------------------
pl.figure()
pl.loglog((1e9)*x_065_par, (1e21)*y_065_par,'b-', label = '[ 6,5]')#A_py_par)
pl.loglog((1e9)*x_090_par, (1e21)*y_090_par,'g-', label = '[ 9,0]')#A_py_par)
pl.loglog((1e9)*x_091_par, (1e21)*y_091_par,'r-', label = '[ 9,1]')#A_py_par)
pl.loglog((1e9)*x_093_par, (1e21)*y_093_par,'c-', label = '[ 9,3]')#A_py_par)
pl.loglog((1e9)*x_290_par, (1e21)*y_290_par,'m-', label = '[29,0]')#A_py_par)
#pl.loglog(x90_A0, y90_A0,'k-' , label = r'GH $\mathcal{A^{(0)}}(\ell)$')
#pl.loglog(x90_A2, y90_A2,'k--', label = r'GH $\mathcal{A^{(2)}}(\ell)$')
pl.xlabel(x_ax)
pl.ylabel(y_ax_par)
pl.title(title('all','all','parallel'))
pl.legend(loc = 'best')
pl.axis([1,500,1e-39,1e4])
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig(svfig('all','all','parallel'))
pl.show()

#--- 65w65 ------------------------------------------------------------------
pl.figure()
pl.loglog((1e9)*x_065_par, (1e21)*y_065_par,'b-', label = '[ 6,5]')#A_py_par)
pl.loglog((1e9)*x_090_par, (1e21)*y_090_par,'g-', label = '[ 9,0]')#A_py_par)
pl.loglog((1e9)*x_091_par, (1e21)*y_091_par,'r-', label = '[ 9,1]')#A_py_par)
pl.loglog((1e9)*x_093_par, (1e21)*y_093_par,'c-', label = '[ 9,3]')#A_py_par)
pl.loglog((1e9)*x_290_par, (1e21)*y_290_par,'m-', label = '[29,0]')#A_py_par)
#pl.loglog(x90_A0, y90_A0,'k-' , label = r'GH $\mathcal{A^{(0)}}(\ell)$')
#pl.loglog(x90_A2, y90_A2,'k--', label = r'GH $\mathcal{A^{(2)}}(\ell)$')
pl.xlabel(x_ax)
pl.ylabel(y_ax_par)
pl.title(title('all','all','parallel'))
pl.legend(loc = 'best')
pl.axis([1,100,1e-4,1e4])
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig(svfig('zoomed','all','parallel'))
pl.show()

#--- 65w65 ------------------------------------------------------------------
pl.figure()
pl.plot((1e9)*x_065_par, abs((1e21)*y_065_par / (1e21)*y_065_par[0]),'b-', label = '[ 6,5]')#A_py_par)
pl.plot((1e9)*x_090_par, abs((1e21)*y_090_par / (1e21)*y_090_par[0]),'g-', label = '[ 9,0]')#A_py_par)
pl.plot((1e9)*x_091_par, abs((1e21)*y_091_par / (1e21)*y_091_par[0]),'r-', label = '[ 9,1]')#A_py_par)
pl.plot((1e9)*x_093_par, abs((1e21)*y_093_par / (1e21)*y_093_par[0]),'c-', label = '[ 9,3]')#A_py_par)
pl.plot((1e9)*x_290_par, abs((1e21)*y_290_par / (1e21)*y_290_par[0]),'m-', label = '[29,0]')#A_py_par)
#pl.loglog(x90_A0, y90_A0,'k-' , label = r'GH $\mathcal{A^{(0)}}(\ell)$')
#pl.loglog(x90_A2, y90_A2,'k--', label = r'GH $\mathcal{A^{(2)}}(\ell)$')
pl.xlabel(x_ax)
pl.ylabel(y_ax_par)
pl.title(title('all','normed','parallel'))
pl.legend(loc = 'best')
#pl.axis([10,1010,14,155])
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig(svfig('all','normed','parallel'))
pl.show()

plsvfig(x_065_par, y_065_par, A_py_par, A_GH_par, x_ax, y_ax_par, title('6' ,'5','parallel'), svfig( '65','65','parallel'))
plsvfig(x_090_par, y_090_par, A_py_par, A_GH_par, x_ax, y_ax_par, title('9' ,'0','parallel'), svfig( '90','90','parallel'))
plsvfig(x_091_par, y_091_par, A_py_par, A_GH_par, x_ax, y_ax_par, title('9' ,'1','parallel'), svfig( '91','91','parallel'))
plsvfig(x_093_par, y_093_par, A_py_par, A_GH_par, x_ax, y_ax_par, title('9' ,'3','parallel'), svfig( '93','93','parallel'))
plsvfig(x_290_par, y_290_par, A_py_par, A_GH_par, x_ax, y_ax_par, title('29','0','parallel'), svfig('290','290','parallel'))

##--- 93w93 ------------------------------------------------------------------
#pl.figure()
#pl.loglog((1e9)*x_093_par, (1e21)*y_093_par,'b-', label = A_py_par)
##pl.loglog(x90_A0, y90_A0,'k-' , label = r'GH $\mathcal{A^{(0)}}(\ell)$')
##pl.loglog(x90_A2, y90_A2,'k--', label = r'GH $\mathcal{A^{(2)}}(\ell)$')
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_par)
#pl.title(title(93,93,'parallel'))
#pl.legend(loc = 'best')
##pl.axis([1e-9,1e-6,1e-24,1e-19])
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.savefig(svfig('93','93','parallel'))
#pl.show()
