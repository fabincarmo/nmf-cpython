Non-Negative Matrix Factorization (NMF) Implementation in C/Python using CBLAS/ATLAS Library
============================================================================================

Cost Function: Kullback-Leibler Divergence.

![equation](http://latex.codecogs.com/gif.latex?D_%7BKL%7D%28X%20%7C%7C%20W%20H%29%20%3D%20%5Csum_%7Bf%7D%20%5Csum_%7Bn%7D%20d%28%5BV%5D_%7Bfn%7D%20%7C%20%5BW%20H%5D_%7Bfn%7D%29)
![equation](http://latex.codecogs.com/gif.latex?d_%7BKL%7D%28x%7Cy%29%20%3D%20x%20%5Clog%5Cleft%28%5Cfrac%7Bx%7D%7By%7D%5Cright%29%20-x%20-%20y)

Requires:

* python-dev or python3-dev
* libatlas3-base 

Building
--------

```
$ make py2 
```
or
```
$ make py3
```

Example
-------
```
$ python run.py
V = 
[[ 1.  1.  1.  1.  1.]
 [ 0.  1.  0.  1.  0.]
 [ 0.  1.  0.  1.  0.]]
W . H ~
[[ 1.  1.  1.  1.  1.]
 [ 0.  1.  0.  1.  0.]
 [ 0.  1.  0.  1.  0.]]
Part 1 =
[[ 0.     0.769  0.     0.769  0.   ]
 [ 0.     1.     0.     1.     0.   ]
 [ 0.     1.     0.     1.     0.   ]]
Part 2 =
[[ 1.     0.231  1.     0.231  1.   ]
 [ 0.     0.     0.     0.     0.   ]
 [ 0.     0.     0.     0.     0.   ]]
```