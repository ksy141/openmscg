#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "force.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() Force *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(compute_pair)
{
    PairList *plist;
    PyArrayObject *dU, *f;
    PyArg_ParseTuple(args, "LOO", &plist, &dU, &f);
    Force::compute_pair(plist, (float*)PyArray_DATA(dU), (float*)PyArray_DATA(f));
    Py_RETURN_NONE;
}

PYAPI(compute_bond)
{
    BondList *blist;
    PyArrayObject *dU, *f;
    PyArg_ParseTuple(args, "LOO", &blist, &dU, &f);
    Force::compute_bond(blist, (float*)PyArray_DATA(dU), (float*)PyArray_DATA(f));
    Py_RETURN_NONE;
}

PYAPI(compute_angle)
{
    BondList *blist;
    PyArrayObject *dU, *f;
    PyArg_ParseTuple(args, "LOO", &blist, &dU, &f);
    Force::compute_angle(blist, (float*)PyArray_DATA(dU), (float*)PyArray_DATA(f));
    Py_RETURN_NONE;
}

PYAPI(compute_dihedral)
{
    BondList *blist;
    PyArrayObject *dU, *f;
    PyArg_ParseTuple(args, "LOO", &blist, &dU, &f);
    Force::compute_dihedral(blist, (float*)PyArray_DATA(dU), (float*)PyArray_DATA(f));
    Py_RETURN_NONE;
}

static PyMethodDef cModPyMethods[] =
{
    {"compute_pair",     compute_pair,     METH_VARARGS, "Compute force"},
    {"compute_bond",     compute_bond,     METH_VARARGS, "Compute force"},
    {"compute_angle",    compute_angle,    METH_VARARGS, "Compute force"},
    {"compute_dihedral", compute_dihedral, METH_VARARGS, "Compute force"},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_force", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_force(void)
{
    import_array();
    return PyModule_Create(&cModPy);
}
