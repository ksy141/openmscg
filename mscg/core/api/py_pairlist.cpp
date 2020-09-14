#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "pair_list.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() PairList *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    float cut, binsize;
    long maxpairs;
    
    PyArg_ParseTuple(args, "ffL", &cut, &binsize, &maxpairs);
    PairList *p = new PairList(cut, binsize, maxpairs);
    return Py_BuildValue("L", p);
}

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

PYAPI(init)
{
    PairList *pair;
    PyArrayObject *types, *exmap;
    PyArg_ParseTuple(args, "LOO", &pair, &types, &exmap);
    
    int *t = (int*)PyArray_DATA(types);
    int natoms = PyArray_DIMS(types)[0];
    
    if(Py_None != (PyObject *)exmap)
        pair->init(t, natoms, (int*)PyArray_DATA(exmap), PyArray_DIMS(exmap)[1]);
    else pair->init(t, natoms, 0, 0);
    
    Py_RETURN_NONE;
}

PYAPI(setup_bins)
{
    PairList *pair;
    PyArrayObject *box;
    
    PyArg_ParseTuple(args, "LO", &pair, &box);
    pair->setup_bins((float*)PyArray_DATA(box));
    Py_RETURN_NONE;
}

PYAPI(build)
{
    PairList *pair;
    PyArrayObject *x;
    
    PyArg_ParseTuple(args, "LO", &pair, &x);
    pair->build((vec3f*)PyArray_DATA(x));
    return Py_BuildValue("L", pair->npairs);
}

PYAPI(fill_page)
{
    PairList *pair;
    int type_id, inext, page_size;
    PyArrayObject *npIndex, *npVector, *npScalar;
    
    PyArg_ParseTuple(args, "LiiiOOO", &pair, &type_id, &inext, &page_size, 
                     &npIndex, &npVector, &npScalar);
    
    int nfill = 0, npairs = pair->npairs;
    
    while(inext<npairs && nfill<page_size)
    {
        while(inext<npairs && pair->tlist[inext]!=type_id) inext++;
        if(inext >= npairs) break;
        
        if(Py_None != (PyObject *)npIndex)
        {
            long *d = (long*)PyArray_DATA(npIndex);
            d[nfill] = pair->ilist[inext];
            d[nfill+page_size] = pair->jlist[inext];
        }
        
        if(Py_None != (PyObject *)npVector)
        {
            double *d = (double*)PyArray_DATA(npVector);
            d[nfill] = pair->dxlist[inext];
            d[nfill+page_size] = pair->dylist[inext];
            d[nfill+page_size+page_size] = pair->dzlist[inext];
        }
        
        if(Py_None != (PyObject *)npScalar)
        {
            double *d = (double*)PyArray_DATA(npScalar);
            d[nfill] = pair->drlist[inext];
        }
        
        nfill++;
        inext++;
    }
    
    return Py_BuildValue("ii", inext, nfill);
}

PYAPI(update_types)
{
    PairList *pair;
    PyArrayObject *npTypes;
    PyArg_ParseTuple(args, "LO", &pair, &npTypes);
    pair->update_types((int*)PyArray_DATA(npTypes));
    Py_RETURN_NONE;
}

PYAPI(get_tid)
{
    int i, j;
    PyArg_ParseTuple(args, "ii", &i, &j);
    return Py_BuildValue("i", pair_tid(i, j));
}

static PyMethodDef cModPyMethods[] =
{
    {"create",       create,       METH_VARARGS, "Create pair-list."},
    {"destroy",      destroy,      METH_VARARGS, "Destroy pair-list."},
    {"init",         init,         METH_VARARGS, "Initialize pair-list."},
    {"setup_bins",   setup_bins,   METH_VARARGS, "Setup verlet-list bins."},
    {"build",        build,        METH_VARARGS, "Build from frame data."},
    {"fill_page",    fill_page,    METH_VARARGS, "Fill a page of pairs."},
    {"get_tid",      get_tid,      METH_VARARGS, "Get type ID of a pair."},
    {"update_types", update_types, METH_VARARGS, "Update pair types."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_pairlist", /* name of module */
    "",         /* module documentation, may be NULL */
    -1,         /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_pairlist(void)
{
    import_array();
    return PyModule_Create(&cModPy);
}
