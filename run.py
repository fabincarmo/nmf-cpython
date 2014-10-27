#!/usr/bin/python

import numpy as np
import nmf
np.set_printoptions(precision=5, suppress=True)

itmax = 500

V = np.array([[1.,1.,1.,1.,1.],[0.,1.,0.,1.,0.],[0.,1.,0.,1.,0.]])
W = np.random.rand(3,2)
H = np.random.rand(2,5)

W2,H2 = nmf.nmf_euc(V,W,H,itmax)

print "V = "
print V
print "W . H ~"
print np.dot(W2,H2)
print np.array(W2)
print np.array(H2)

W2,H2 = nmf.nmf_euc_sparse(V,W,H,0.2,itmax)

print "V = "
print V
print "W . H ~"
print np.dot(W2,H2)
print np.array(W2)
print np.array(H2)

