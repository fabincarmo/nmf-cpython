#include <Python.h>
#include <cblas.h>

#define MULT(A,B,C,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){C[i*m+j] = A[i*m+j] * B[i*m+j];}}
#define DIV(A,B,C,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){C[i*m+j] = A[i*m+j] / B[i*m+j];}}
#define MULTM(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, r, 1.0, A, r, B, m, 0.0, C, m);
#define MULTMtnt(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, r, 1.0, A, n, B, m, 0.0, C, m);
#define MULTMntt(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, r, 1.0, A, r, B, r, 0.0, C, m);

static PyObject* double_nmf(PyObject *self, PyObject *args);

static PyMethodDef nmfMethods[] = {
    {"nmf",
        double_nmf,
        METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

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

static double * NormaColuna(double * M, const int n, const int m){
	int i, j;
	double norma;
	for(i=0;i<m;i++){
		norma = cblas_dnrm2(n,&M[i],m);
		for(j=0;j<n;j++){
			M[i+j*m] = M[i+j*m] / norma;
		}
	}
	return M;
}

static PyObject* double_nmf(PyObject *self, PyObject *args)
{
	int it, n, r, m, i ,j;
	const int itmax;
	double *V, *W, *H, *R, *aux1, *aux2, *Wn;
	PyObject *t1, *t2, *t3;

    PyArg_ParseTuple(args, "OOOi", &t1, &t2, &t3, &itmax);

	n = PySequence_Size(t1);
	m = PySequence_Size(PySequence_GetItem(t1,0));
	r = PySequence_Size(t3);

	V = list2v(t1);
	W = list2v(t2);
	H = list2v(t3);

	R = malloc(sizeof(double) * n*m);
	Wn = malloc(sizeof(double) * n*r);

	for(it=0; it<itmax; it++){
		Wn = NormaColuna(W,n,r);
		MULTM(W, H, R, n, r, m);
		aux1 = malloc(sizeof(double) * r*m);
		aux2 = malloc(sizeof(double) * r*m);
		MULTMtnt(Wn, V, aux1, r,n,m);
		MULTMtnt(Wn, R, aux2, r,n,m);
		DIV(aux1, aux2, aux1, r,m);
		MULT(H, aux1, H, r,m);
		free(aux1);
		free(aux2);

		aux1 = malloc(sizeof(double) * n*r);
		aux2 = malloc(sizeof(double) * n*r);
		MULTM(W, H, R, n,r,m);
		MULTMntt(V, H, aux1, n,m,r);
		MULTMntt(R, H, aux2, n,m,r);
		DIV(aux1, aux2, aux1, n,r);
		MULT(W, aux1, W, n,r);
		free(aux1);
		free(aux2);

	}

	t1 = v2list(V,n,m);
	t2 = v2list(W,n,r);
	t3 = v2list(H,r,m);

    return Py_BuildValue("OO", t2, t3);
}

PyMODINIT_FUNC initnmf(void)
{
    PyObject *m;

    m = Py_InitModule("nmf", nmfMethods);
    if (m == NULL) {
        return;
    }
}
