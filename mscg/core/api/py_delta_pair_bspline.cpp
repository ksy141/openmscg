#include <Python.h>
#include "delta_pair_bspline.h"
#include "pair_list.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() DeltaPairBSpline *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    PairList *pair;
    int tid, order;
    double res, xmin; 
    PyArg_ParseTuple(args, "Liidd", &pair, &tid, &order, &res, &xmin);
    printf("DeltaPairBSpline %d %lf %lf\n", order, res, xmin);
    DeltaPairBSpline *p = new DeltaPairBSpline(pair, tid, order, res, xmin);
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
    DeltaPairBSpline *p;
    double ddx_factor;
    PyArg_ParseTuple(args, "Ld", &p, &ddx_factor);
    p->setup_cache(ddx_factor);
    Py_RETURN_NONE;
}

PYAPI(compute)
{
    GETPTR();
    p->compute();
    Py_RETURN_NONE;
}

PYAPI(get_spline)
{
    GETPTR();
    return Py_BuildValue("iiifff", p->order, p->nbreak, p->ncoeff, p->xmin, p->xmax, p->resolution);
}

static PyMethodDef cModPyMethods[] =
{
    {"create",      create,      METH_VARARGS, "Create delta object."},
    {"destroy",     destroy,     METH_VARARGS, "Destroy delta object."},
    {"setup_cache", setup_cache, METH_VARARGS, "Setup cache acceleration."},
    {"compute",     compute,     METH_VARARGS, "Compute coefficients."},
    {"get_spline",  get_spline,  METH_VARARGS, "Get spline parameters."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_delta_pair_bspline", /* name of module */
    "",         /* module documentation, may be NULL */
    -1,         /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_delta_pair_bspline(void)
{
    return PyModule_Create(&cModPy);
}
