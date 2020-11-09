#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "bond_list.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() BondList *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    PyArrayObject *npBT, *npAT, *npDT, *npBonds, *npAngles, *npDihedrals;
    PyArg_ParseTuple(args, "OOOOOO", &npBT, &npBonds, &npAT, &npAngles, &npDT, &npDihedrals);
    
    int nbonds = 0, nangles = 0, ndihedrals = 0;
    int *bt = 0, *at = 0, *dt = 0;
    vec2i *bonds = NULL;
    vec3i *angles = NULL;
    vec4i *dihedrals = NULL;
    
    if(Py_None != (PyObject *)npBonds)
    {
        nbonds = PyArray_DIMS(npBonds)[0];
        bt = (int*)PyArray_DATA(npBT);
        bonds = (vec2i*)PyArray_DATA(npBonds);
    }
    
    if(Py_None != (PyObject *)npAngles)
    {
        nangles = PyArray_DIMS(npAngles)[0];
        at = (int*)PyArray_DATA(npAT);
        angles = (vec3i*)PyArray_DATA(npAngles);
    }
    
    if(Py_None != (PyObject *)npDihedrals)
    {
        ndihedrals = PyArray_DIMS(npDihedrals)[0];
        dt = (int*)PyArray_DATA(npDT);
        dihedrals = (vec4i*)PyArray_DATA(npDihedrals);
    }
    
    BondList *p = new BondList(nbonds, bt, bonds, nangles, at, angles, ndihedrals, dt, dihedrals);
    return Py_BuildValue("L", p);
}

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

PYAPI(build)
{
    BondList *p;
    PyArrayObject *npBox, *npX;
    
    PyArg_ParseTuple(args, "LOO", &p, &npBox, &npX);
    p->build((float*)PyArray_DATA(npBox), (vec3f*)PyArray_DATA(npX));
    Py_RETURN_NONE;
}

PYAPI(get_scalar)
{
    BondList *p;
    int target;
    PyArrayObject *npData;
    PyArg_ParseTuple(args, "LiO", &p, &target, &npData);
    
    int n = (target==0?p->nbonds:(target==1?p->nangles:p->ndihedrals));
    float *pd = (target==0?p->dr_bond:(target==1?p->theta_angle:p->phi_dihedral));
    
    float *des = (float*)PyArray_DATA(npData);
    for(int i=0; i<n; i++) des[i] = pd[i];
    Py_RETURN_NONE;
}

static PyMethodDef cModPyMethods[] =
{
    {"create",     create,     METH_VARARGS, "Create pair-list."},
    {"destroy",    destroy,    METH_VARARGS, "Destroy pair-list."},
    {"build",      build,      METH_VARARGS, "Build from frame data."},
    {"get_scalar", get_scalar, METH_VARARGS, "Get bonding values."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_bondlist", /* name of module */
    "",         /* module documentation, may be NULL */
    -1,         /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_bondlist(void)
{
    import_array();
    return PyModule_Create(&cModPy);
}
