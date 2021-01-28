#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "tables.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() Tables *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    int maxtable;
    PyArg_ParseTuple(args, "i", &maxtable);
    Tables *p = new Tables(maxtable);
    return Py_BuildValue("L", p);
}

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

void convert_table(PyObject* p, Tables::Table *tb)
{
    tb->min = (float)PyFloat_AsDouble(PyDict_GetItemString(p, "min"));
    tb->inc = (float)PyFloat_AsDouble(PyDict_GetItemString(p, "inc"));
    PyArrayObject *np_table = (PyArrayObject*)PyDict_GetItemString(p, "efac");
    tb->n = PyArray_DIMS(np_table)[0];
    tb->efac = (float*)PyArray_DATA(np_table);
    
    np_table = (PyArrayObject*)PyDict_GetItemString(p, "ffac");
    tb->ffac = (float*)PyArray_DATA(np_table);
}

PYAPI(set_table)
{
    Tables *p;
    int tid;
    PyObject *tbl;
    PyArg_ParseTuple(args, "LiO", &p, &tid, &tbl);
    convert_table(tbl, p->tbls + tid);
    Py_RETURN_NONE;
}

PYAPI(compute)
{
    Tables *p;
    long n;
    PyArrayObject *types, *scalars, *U, *dU;
    PyArg_ParseTuple(args, "LLOOOO", &p, &n, &types, &scalars, &U, &dU);
    
    p->compute(n,
        (int*)PyArray_DATA(types),
        (float*)PyArray_DATA(scalars),
        (float*)PyArray_DATA(U),
        (float*)PyArray_DATA(dU));
    
    Py_RETURN_NONE;
}
/*
PYAPI(compute_bond)
{
    int type_id;
    Force *force;
    BondList *blist;
    PyArrayObject *f;
    PyObject *tbl;
    PyArg_ParseTuple(args, "LiLOO", &force, &type_id, &blist, &f, &tbl);
    
    ForceTable tb;
    convert_table(tbl, &tb);
    force->compute_bond(type_id, blist, (float*)PyArray_DATA(f), &tb);
    Py_RETURN_NONE;
}

PYAPI(compute_angle)
{
    int type_id;
    Force *force;
    BondList *blist;
    PyArrayObject *f;
    PyObject *tbl;
    PyArg_ParseTuple(args, "LiLOO", &force, &type_id, &blist, &f, &tbl);
    
    ForceTable tb;
    convert_table(tbl, &tb);
    force->compute_angle(type_id, blist, (float*)PyArray_DATA(f), &tb);
    Py_RETURN_NONE;
}

PYAPI(compute_dihedral)
{
    int type_id;
    Force *force;
    BondList *blist;
    PyArrayObject *f;
    PyObject *tbl;
    PyArg_ParseTuple(args, "LiLOO", &force, &type_id, &blist, &f, &tbl);
    
    ForceTable tb;
    convert_table(tbl, &tb);
    force->compute_dihedral(type_id, blist, (float*)PyArray_DATA(f), &tb);
    Py_RETURN_NONE;
}
*/
static PyMethodDef cModPyMethods[] =
{
    {"create",    create,    METH_VARARGS, "Create tables module."},
    {"destroy",   destroy,   METH_VARARGS, "Destroy tables module."},
    {"set_table", set_table, METH_VARARGS, "Setup a table."},
    {"compute",   compute,   METH_VARARGS, "Compute dU values."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_tables", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_tables(void)
{
    import_array();
    return PyModule_Create(&cModPy);
}
