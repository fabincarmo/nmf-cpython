#!/usr/bin/python

import random
import numpy as np
import nmf
from time import time

n = 512
r = 64
m = 100
itmax = 150000
V = np.array([[1,1,1,1,1],[0,1,0,1,0],[0,1,0,1,0]])
W = np.array([[.5,.2],[.1,.5],[.1,.5]])
H = np.array([[.5,.1,.1,.5,.5],[.1,.5,.8,.9,.5]])

tempos=[]
for i in range(1):
    # V = np.random.rand(n,m)
    # W = np.random.rand(n,r)
    # H = np.random.rand(r,m)
    t1 = time()
    W2 = W
    H2 = H
    for j in range(itmax):
        Wn = W2 / np.sum(W2**2, axis=0)**(0.5)
        H2 = H2 * (np.dot(Wn.T, V)) / (np.dot(Wn.T, np.dot(Wn,H2)))
        W2 = W2 * (np.dot(V, H2.T)) / (np.dot(np.dot(W2,H2),H2.T))
    t2 = time() - t1
    print t2,
    tempos.append(t2)
print np.mean(tempos)
print np.dot(W2,H2)

tempos=[]
for i in range(1):
    # V = np.random.rand(n,m)
    # W = np.random.rand(n,r)
    # H = np.random.rand(r,m)
    t1 = time()
    W2,H2 = nmf.nmf(V,W,H,itmax)
    t2 = time() - t1
    print t2,
    tempos.append(t2)
print np.mean(tempos)
print np.dot(W2,H2)
