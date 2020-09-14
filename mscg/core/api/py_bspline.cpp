#include <Python.h>
#include "bspline.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() BSpline *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    int order;
    double resolution, xmin, xmax;
    
    PyArg_ParseTuple(args, "iddd", &order, &resolution, &xmin, &xmax);
    BSpline *p = new BSpline(order, resolution, xmin, xmax);
    return Py_BuildValue("L", p);
}

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

PYAPI(interp)
{
    BSpline *p;
    double xmin, dx;
    int n;
    PyObject *pList;
    
    PyArg_ParseTuple(args, "LddiO", &p, &xmin, &dx, &n, &pList);
    
    double *input = new double[p->ncoeff];
    double *output = new double[n];
    
    for(int i=0; i<p->ncoeff; i++) 
    {
        PyObject *pItem = PyList_GetItem(pList, i);
        input[i] = PyFloat_AsDouble(pItem);
    }
    
    p->eval(input, output, xmin, dx, n);
    
    PyObject *z = PyList_New(n);
    for(int i=0; i<n; i++) PyList_SetItem(z, i, Py_BuildValue("f", output[i]));
    
    delete [] input;
    delete [] output;
    return z;
}

static PyMethodDef cModPyMethods[] =
{
    {"create",  create,  METH_VARARGS, "Create B-spline object."},
    {"destroy", destroy, METH_VARARGS, "Destroy B-spline object."},
    {"interp",  interp,  METH_VARARGS, "Interpolate spline."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_bspline", /* name of module */
    "",         /* module documentation, may be NULL */
    -1,         /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_bspline(void)
{
    return PyModule_Create(&cModPy);
}
