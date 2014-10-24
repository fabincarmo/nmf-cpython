#!/usr/bin/python

import numpy as np
import nmf
np.set_printoptions(precision=5, suppress=True)

itmax = 1500

V = np.array([[1.,1.,1.,1.,1.],[0.,1.,0.,1.,0.],[0.,1.,0.,1.,0.]])
W = np.random.rand(3,2)
H = np.random.rand(2,5)

W2,H2 = nmf.nmf(V,W,H,itmax)

print "V = "
print V
print "W . H ~"
print np.dot(W2,H2)
