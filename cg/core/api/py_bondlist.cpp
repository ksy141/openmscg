#include <Python.h>
#include "bond_list.h"
#include "traj.h"
#include "topology.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() BondList *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    Topology *top;
    
    PyArg_ParseTuple(args, "L", &top);
    BondList *p = new BondList(top);
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
    Traj *traj;
    
    PyArg_ParseTuple(args, "LL", &p, &traj);
    p->build(traj);
    Py_RETURN_NONE;
}

PYAPI(get_bonds)
{
    GETPTR();
    Topology *top = p->top;
    int count = top->nbonds;
    
    PyObject *tlist = PyList_New(count);
    PyObject *ilist = PyList_New(count);
    PyObject *jlist = PyList_New(count);
    PyObject *rlist = PyList_New(count);
    
    for(int i=0; i<count; i++)
    {
        PyList_SetItem(tlist, i, Py_BuildValue("i", top->bond_types[i]));
        PyList_SetItem(ilist, i, Py_BuildValue("i", top->bond_atom1[i]));
        PyList_SetItem(jlist, i, Py_BuildValue("i", top->bond_atom2[i]));
        PyList_SetItem(rlist, i, Py_BuildValue("f", p->dr_bond[i]));
    }
    
    PyObject *z = Py_BuildValue("(O,O,O,O)", tlist, ilist, jlist, rlist);
    Py_XDECREF(tlist); Py_XDECREF(ilist); Py_XDECREF(jlist); Py_XDECREF(rlist);
    return z;
}

PYAPI(get_angles)
{
    GETPTR();
    Topology *top = p->top;
    int count = top->nangls;
    
    PyObject *tlist = PyList_New(count);
    PyObject *ilist = PyList_New(count);
    PyObject *jlist = PyList_New(count);
    PyObject *klist = PyList_New(count);
    PyObject *rlist = PyList_New(count);
    
    for(int i=0; i<count; i++)
    {
        PyList_SetItem(tlist, i, Py_BuildValue("i", top->angl_types[i]));
        PyList_SetItem(ilist, i, Py_BuildValue("i", top->angl_atom1[i]));
        PyList_SetItem(jlist, i, Py_BuildValue("i", top->angl_atom2[i]));
        PyList_SetItem(klist, i, Py_BuildValue("i", top->angl_atom3[i]));
        PyList_SetItem(rlist, i, Py_BuildValue("f", p->theta_angl[i]));
    }
    
    PyObject *z = Py_BuildValue("(O,O,O,O,O)", tlist, ilist, jlist, klist, rlist);
    Py_XDECREF(tlist); Py_XDECREF(ilist); Py_XDECREF(jlist); Py_XDECREF(klist); Py_XDECREF(rlist);
    return z;
}

static PyMethodDef cModPyMethods[] =
{
    {"create",     create,     METH_VARARGS, "Create pair-list."},
    {"destroy",    destroy,    METH_VARARGS, "Destroy pair-list."},
    {"build",      build,      METH_VARARGS, "Build from frame data."},
    {"get_bonds",  get_bonds,  METH_VARARGS, "Get bonds data."},
    {"get_angles", get_angles, METH_VARARGS, "Get angles data."},
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
    return PyModule_Create(&cModPy);
}
