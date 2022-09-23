#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "nb3b.h"
#include "pair_list.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() NB3B *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

PYAPI(create_sw)
{
    int tid_i, tid_ij, tid_ik;
    double gamma_ij, a_ij, gamma_ik, a_ik;
    double lambda, cos0;

    PyArg_ParseTuple(args, "iiidddddd", &tid_i, &tid_ij, &tid_ik,
        &lambda, &cos0, &gamma_ij, &a_ij, &gamma_ik, &a_ik);

    NB3B_SW *p = new NB3B_SW(tid_i, tid_ij, tid_ik,
        lambda, cos0, gamma_ij, a_ij, gamma_ik, a_ik);

    return Py_BuildValue("L", p);
}

PYAPI(compute_sw)
{
    NB3B_SW *p;
    PairList *plist;
    PyArrayObject *U, *dU, *F;

    PyArg_ParseTuple(args, "LLOOO", &p, &plist, &U, &dU, &F);

    p->compute(plist, (float*)PyArray_DATA(U),
        (float*)PyArray_DATA(dU), (float*)PyArray_DATA(F));

    Py_RETURN_NONE;
}

static PyMethodDef cModPyMethods[] =
{
    {"destroy",    destroy,    METH_VARARGS, "Destroy NB3B module."},
    {"create_sw",  create_sw,  METH_VARARGS, "Create SW."},
    {"compute_sw", compute_sw, METH_VARARGS, "Compute dU values with SW."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_nb3b", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_nb3b(void)
{
    import_array();
    return PyModule_Create(&cModPy);
}
