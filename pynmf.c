#include "array_list.c"
#include "nmf.c"

static PyObject* double_nmf_kl(PyObject *self, PyObject *args);
static PyObject* double_nmf_kl_h(PyObject *self, PyObject *args);

static PyMethodDef nmf_Methods[] = {
    {"nmf_kl", double_nmf_kl, METH_VARARGS, ""},
	{"nmf_kl_h", double_nmf_kl_h, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};


static PyObject* double_nmf_kl(PyObject *self, PyObject *args)
{
	int n, m, r;
	const int itmax;
	double *V, *W, *H;
	PyObject *t1, *t2, *t3;

    PyArg_ParseTuple(args, "OOOi", &t1, &t2, &t3, &itmax);


	V = list2v(t1);
	W = list2v(t2);
	H = list2v(t3);

	n = PySequence_Size(t1);
	m = PySequence_Size(PySequence_GetItem(t1,0));
	r = PySequence_Size(t3);

	nmf_kl(V, W, H, n, r, m, itmax);

	t2 = v2list(W,n,r);
	t3 = v2list(H,r,m);

    return Py_BuildValue("OO", t2, t3);
}

static PyObject* double_nmf_kl_h(PyObject *self, PyObject *args)
{
	int n, r, m;
	const int itmax;
	double *V, *W, *H;
	PyObject *t1, *t2, *t3;

    PyArg_ParseTuple(args, "OOOi", &t1, &t2, &t3, &itmax);

	n = PySequence_Size(t1);
	m = PySequence_Size(PySequence_GetItem(t1,0));
	r = PySequence_Size(t3);

	V = list2v(t1);
	W = list2v(t2);
	H = list2v(t3);

	nmf_kl_h(V, W, H, n, r, m, itmax);

	t2 = v2list(W,n,r);
	t3 = v2list(H,r,m);

    return Py_BuildValue("OO", t2, t3);
}

#if PY_VERSION_HEX > 0x03000000

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "nmf",
    NULL,
    -1,
    nmf_Methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_nmf(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

#else

PyMODINIT_FUNC initnmf(void)
{
    PyObject *m;

    m = Py_InitModule("nmf", nmf_Methods);
    if (m == NULL) {
        return;
    }
}

#endif
