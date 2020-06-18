#include <Python.h>
#include <numpy/arrayobject.h>

#include "pair_list.h"
#include "traj.h"
#include "topology.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() PairList *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    Topology *top;
    int natoms;
    
    PyArg_ParseTuple(args, "Li", &top, &natoms);
    PairList *p = new PairList(top, natoms);
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
    float cut, binsize;
    PyArg_ParseTuple(args, "Lff", &pair, &cut, &binsize);
    pair->init(cut, binsize);
    Py_RETURN_NONE;
}

PYAPI(setup_bins)
{
    PairList *pair;
    Traj *traj;
    
    PyArg_ParseTuple(args, "LL", &pair, &traj);
    pair->setup_bins(traj);
    Py_RETURN_NONE;
}

PYAPI(build)
{
    PairList *pair;
    Traj *traj;
    int reset_bins;
    
    PyArg_ParseTuple(args, "LLp", &pair, &traj, &reset_bins);
    pair->build(traj, reset_bins == 1);
    return Py_BuildValue("L", pair->npairs);
}

PYAPI(get_pairs)
{
    PairList *pair;
    int start, count;
    
    PyArg_ParseTuple(args, "Lii", &pair, &start, &count);
    if(start + count > pair->npairs) count = pair->npairs - start;
    if(count<=0) return Py_BuildValue("[[],[],[],[]]");
    
    PyObject *tlist = PyList_New(count);
    PyObject *ilist = PyList_New(count);
    PyObject *jlist = PyList_New(count);
    PyObject *rlist = PyList_New(count);
    
    for(int i=0; i<count; i++)
    {
        PyList_SetItem(tlist, i, Py_BuildValue("i", pair->tlist[i+start]));
        PyList_SetItem(ilist, i, Py_BuildValue("i", pair->ilist[i+start]));
        PyList_SetItem(jlist, i, Py_BuildValue("i", pair->jlist[i+start]));
        PyList_SetItem(rlist, i, Py_BuildValue("f", pair->drlist[i+start]));
    }
    
    return Py_BuildValue("(N,N,N,N)", tlist, ilist, jlist, rlist);
}

PYAPI(fill_page)
{
    PairList *pair;
    int inext, page_size, type_id;
    PyArrayObject *npIndex, *npVector, *npScalar;
    
    PyArg_ParseTuple(args, "LiiiOOO", &pair, &type_id, &inext, &page_size, &npIndex, &npVector, &npScalar);
    int nfill = 0, npairs = pair->npairs;
    
    while(inext<npairs && nfill<page_size)
    {
        while(inext<npairs && pair->tlist[inext]!=type_id) inext++;
        if(inext >= npairs) break;
        
        if(Py_None != (PyObject *)npIndex)
        {
            long *d = (long*)(npIndex->data);
            d[nfill] = pair->ilist[inext];
            d[nfill+page_size] = pair->jlist[inext];
        }
        
        if(Py_None != (PyObject *)npVector)
        {
            double *d = (double*)(npVector->data);
            d[nfill] = pair->dxlist[inext];
            d[nfill+page_size] = pair->dylist[inext];
            d[nfill+page_size+page_size] = pair->dzlist[inext];
        }
        
        if(Py_None != (PyObject *)npScalar)
        {
            double *d = (double*)(npScalar->data);
            d[nfill] = pair->drlist[inext];
        }
        
        nfill++;
        inext++;
    }
    
    return Py_BuildValue("ii", inext, nfill);
}

static PyMethodDef cModPyMethods[] =
{
    {"create",     create,     METH_VARARGS, "Create pair-list."},
    {"destroy",    destroy,    METH_VARARGS, "Destroy pair-list."},
    {"init",       init,       METH_VARARGS, "Initialize pair-list."},
    {"setup_bins", setup_bins, METH_VARARGS, "Setup verlet-list bins."},
    {"build",      build,      METH_VARARGS, "Build from frame data."},
    {"get_pairs",  get_pairs,  METH_VARARGS, "Get pairs data."},
    {"fill_page",  fill_page,   METH_VARARGS, "Fill a page of pairs."},
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
