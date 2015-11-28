#!/usr/bin/python

import numpy as np
import nmf
np.set_printoptions(precision=3, suppress=True)

itmax = 100

V = np.array  ([[1., 1., 1., 1., 1.],
                [0., 1., 0., 1., 0.],
                [0., 1., 0., 1., 0.]])
W = np.random.rand(3, 2)
H = np.random.rand(2, 5)

W, H = nmf.nmf_kl(V, W, H, itmax)
W = np.array(W); H = np.array(H)

print("V = ")
print(V)
print("W . H ~")
print(np.dot(W, H))
print("Part 1 =")
print(np.multiply(W[:,0,np.newaxis],H[np.newaxis,0,:]))
print("Part 2 =")
print(np.multiply(W[:,1,np.newaxis],H[np.newaxis,1,:]))
