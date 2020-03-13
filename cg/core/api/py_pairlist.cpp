#include <Python.h>
#include "pair_list.h"
#include "traj.h"
#include "topology.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() PairList *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    Topology *top;
    Traj *traj;
    
    PyArg_ParseTuple(args, "LL", &top, &traj);
    PairList *p = new PairList(top, traj);
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
    GETPTR();
    p->setup_bins();
    Py_RETURN_NONE;
}

PYAPI(build)
{
    PairList *pair;
    int reset_bins;
    PyArg_ParseTuple(args, "Lp", &pair, &reset_bins);
    pair->build(reset_bins == 1);
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
    
    PyObject *z = Py_BuildValue("(O,O,O,O)", tlist, ilist, jlist, rlist);
    Py_XDECREF(tlist); Py_XDECREF(ilist); Py_XDECREF(jlist); Py_XDECREF(rlist);
    return z;
}

static PyMethodDef cModPyMethods[] =
{
    {"create",     create,     METH_VARARGS, "Create pair-list."},
    {"destroy",    destroy,    METH_VARARGS, "Destroy pair-list."},
    {"init",       init,       METH_VARARGS, "Initialize pair-list."},
    {"setup_bins", setup_bins, METH_VARARGS, "Setup verlet-list bins."},
    {"build",      build,      METH_VARARGS, "Build from frame data."},
    {"get_pairs",  get_pairs,  METH_VARARGS, "Get pairs data."},
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
    return PyModule_Create(&cModPy);
}
