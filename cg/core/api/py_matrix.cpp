#include <Python.h>
#include "matrix.h"
#include "table.h"
#include "traj.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() Matrix *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    Matrix *p = new Matrix();
    return Py_BuildValue("L", p);
}

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

PYAPI(add_table)
{
    Matrix *mtx; 
    Table *tbl;
    PyArg_ParseTuple(args, "LL", &mtx, &tbl);
    mtx->add_table(tbl);
    Py_RETURN_NONE;
}

PYAPI(setup)
{
    Matrix *mtx; 
    int natoms;
    PyArg_ParseTuple(args, "Li", &mtx, &natoms);
    mtx->setup(natoms);
    Py_RETURN_NONE;
}

PYAPI(multiplyadd)
{
    Matrix *mtx; 
    Traj *trj;
    PyArg_ParseTuple(args, "LL", &mtx, &trj);
    mtx->multiplyadd((float*)(trj->f));
    Py_RETURN_NONE;
}

PYAPI(reset) { GETPTR(); p->reset(); Py_RETURN_NONE; }
PYAPI(solve) { GETPTR(); p->solve(); Py_RETURN_NONE; }

static PyMethodDef cModPyMethods[] =
{
    {"create",      create,      METH_VARARGS, "Create matrix object."},
    {"destroy",     destroy,     METH_VARARGS, "Destroy matrix object."},
    {"add_table",   add_table,   METH_VARARGS, "Add a table."},
    {"setup",       setup,       METH_VARARGS, "Setup matrix workspace."},
    {"reset",       reset,       METH_VARARGS, "Zero coefficient matrix."},
    {"solve",       solve,       METH_VARARGS, "Solve SVD."},
    {"multiplyadd", multiplyadd, METH_VARARGS, "C = M**T * M + C."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_matrix", /* name of module */
    "",         /* module documentation, may be NULL */
    -1,         /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_matrix(void)
{
    return PyModule_Create(&cModPy);
}


