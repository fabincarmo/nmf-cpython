Non-Negative Matrix Factorization (NMF) Implementation in C/Python using CBLAS/ATLAS Library
============================================================================================

Cost Function: Kullback-Leibler Divergence.
![equation](http://www.sciweavers.org/tex2img.php?eq=D_%7BKL%7D%28X%20%7C%7C%20%5Ctilde%20X%29%20%3D%20%5Csum_%7Bi%2Cj%7D%20X_%7Bi%2Cj%7D%20%5Clog%20%5Cfrac%7BX_%7Bi%2Cj%7D%7D%7B%5Ctilde%20X_%7Bi%2Cj%7D%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

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