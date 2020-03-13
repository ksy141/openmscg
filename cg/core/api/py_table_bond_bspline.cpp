#include <Python.h>
#include "table_bond_bspline.h"
#include "topology.h"
#include "bond_list.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() TableBondBSpline *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    BondList *lst;
    int tid, order;
    double res, xmin, xmax;
    PyArg_ParseTuple(args, "Liiddd", &lst, &tid, &order, &res, &xmin, &xmax);
    TableBondBSpline *p = new TableBondBSpline(lst, tid, order, res, xmin, xmax);
    return Py_BuildValue("L", p);
}

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

PYAPI(setup_cache)
{
    TableBondBSpline *p;
    double dx_factor;
    PyArg_ParseTuple(args, "Ld", &p, &dx_factor);
    p->setup_cache(dx_factor);
    Py_RETURN_NONE;
}

PYAPI(compute)
{
    GETPTR();
    p->compute();
    Py_RETURN_NONE;
}

PYAPI(dump)
{
    TableBondBSpline *p;
    double xmin, dx;
    int n;
    PyArg_ParseTuple(args, "Lddi", &p, &xmin, &dx, &n);

    double *tbl = new double[n];
    p->dump(tbl, xmin, dx, n);

    PyObject *z  = PyList_New(n);
    for(int i=0; i<n; i++) PyList_SetItem(z,  i, Py_BuildValue("d", tbl[i]));
    return z;
}

static PyMethodDef cModPyMethods[] =
{
    {"create",      create,      METH_VARARGS, "Create table object."},
    {"destroy",     destroy,     METH_VARARGS, "Destroy table object."},
    {"setup_cache", setup_cache, METH_VARARGS, "Setup cache acceleration."},
    {"compute",     compute,     METH_VARARGS, "Compute coefficients."},
    {"dump",        dump,        METH_VARARGS, "Dump result table."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_table_bond_bspline", /* name of module */
    "",         /* module documentation, may be NULL */
    -1,         /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_table_bond_bspline(void)
{
    return PyModule_Create(&cModPy);
}
