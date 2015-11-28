Non-Negative Matrix Factorization (NMF) Implementation in C/Python using CBLAS/ATLAS Library
============================================================================================

Cost Function: Kullback-Leibler Divergence.

![equation](http://latex.codecogs.com/gif.latex?D_{KL}(X&space;||&space;\tilde&space;X)&space;=&space;\sum_{i,j}&space;X_{i,j}&space;\log&space;\frac{X_{i,j}}{\tilde&space;X_{i,j}})

Requires:

python-dev or python3-dev
libatlas3-base 

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
```