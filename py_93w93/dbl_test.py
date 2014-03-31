#!/usr/bin/env python

import numpy as np
#import scipy.integrate as integrate
from scipy.integrate import quad

#def integrand(t,n,x,b):
#    return (1./np.sqrt(x*x - 1.)) * t * (1./np.sqrt(t*t + 1.)) * np.exp(-2. * x * np.sqrt(t*t +1.)) 
#
#def expint(n,x,b):
#    return quad(integrand, 1, np.inf, args=(n, x,b))[0]
#
#result = quad(lambda x: expint(n,x,b), 0, np.inf)

#ns = np.arange(1.,5.) 
#bs = np.arange(11.,15.) 
#
#for j,n in enumerate(ns):
#	for k,b in enumerate(bs):
#		result = quad(lambda x: expint(n, x,b), 0, np.inf)
#		print result
#

def integrand(t,x,n,b):
    return (1./np.sqrt(x*x - 1.)) * t * (1./np.sqrt(t*t + 1.)) * np.exp(-n * x * b * np.sqrt(t*t +1.)) 

hs = np.linspace(1.1,11.1,10)
ts = np.linspace(1,11,10)

n = 2.
b = 1.0e-9
f = np.zeros(shape=(len(ts),len(hs)))

for i,time in enumerate(ts):
    for j,h in enumerate(hs):
	f[i,j] = integrand(time,h,2.0,1.0e-9)
F = np.sum(f, axis = 0)
FF = np.sum(F)

print 'double integral  = %2.2',%FF

#def expint(x,n,b):
#    return quad(integrand, 1.1, 10.0, args=(x, n,b))[0]
#
#result = quad(lambda x: expint(x,n,b), 0, 1.0)# np.inf)
#print result
