#!/usr/bin/env python

import numpy as np
#import scipy.integrate as integrate
from scipy.integrate import quad

def integrand(t,n,x,b):
    return b * np.exp(-x*t) / t**n

def expint(n,x,b):
    return quad(integrand, 1, np.inf, args=(n, x,b))[0]


ns = np.arange(1.,5.) 
bs = np.arange(11.,15.) 

for j,n in enumerate(ns):
	for k,b in enumerate(bs):
		result = quad(lambda x: expint(n, x,b), 0, np.inf)
		print result


