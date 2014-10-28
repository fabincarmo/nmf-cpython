#include <Python.h>

/*
 * Converte Python List em C Array
 */
static double* list2v(PyObject *t)
{
	int i, j, n, m;
	double *v;
	PyObject *linha;

	n = PySequence_Size(t);
	m = PySequence_Size(PySequence_GetItem(t,0));

	v = malloc( sizeof(double) * n*m);
	for(i=0; i<n; i++){
		linha = PySequence_GetItem(t,i);
		for(j=0; j<m; j++){
			v[ j+i*m ] = PyFloat_AsDouble(PySequence_GetItem(linha,j));
		}
	}
	return v;
}

/*
 * Converte C Array em Python List
 */
static PyObject* v2list(double *v, const int n, const int m)
{
	int i, j;
	PyObject *aux, *l;

	l = PyList_New(n);
	for(i=0; i<n; i++){
		aux = PyList_New(m);
		for(j=0; j<m; j++){
			PyList_SetItem(aux, j, PyFloat_FromDouble(v[j+i*m]));
		}
		PyList_SetItem(l, i, aux);
	}
	return l;
}

