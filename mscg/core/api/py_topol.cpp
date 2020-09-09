#include <Python.h>
#include "topology.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETPTR() Topology *p; PyArg_ParseTuple(args, "L", &p)

PYAPI(create)
{
    int ntypes;
    PyArg_ParseTuple(args, "i", &ntypes);
    Topology *p = new Topology(ntypes);
    return Py_BuildValue("L", p);
}

PYAPI(destroy)
{
    GETPTR();
    delete p;
    Py_RETURN_NONE;
}

PYAPI(add_atoms)
{
    Topology *top;
    PyObject *pList;
    PyArg_ParseTuple(args, "LO", &top, &pList);
    
    int n = PyList_Size(pList);
    int types[1000];
    int t = 0;
    
    for(int i=0; i<n; i++)
    {
        types[t++] = static_cast<int>(PyLong_AsLong(PyList_GetItem(pList, i)));
        
        if(t==1000)
        {
            top->add_atoms(t, types);
            t = 0;
        }
    }
    
    if(t>0) top->add_atoms(t, types);
    return Py_BuildValue("i", top->natoms);
}

PYAPI(add_bonds)
{
    Topology *top;
    PyObject *pType, *pAtom1, *pAtom2;
    PyArg_ParseTuple(args, "LOOO", &top, &pType, &pAtom1, &pAtom2);
    
    int t = 0, n = PyList_Size(pType);
    int types[1000], atom1[1000], atom2[1000];
    
    for(int i=0; i<n; i++)
    {
        atom1[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom1, i)));
        atom2[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom2, i)));
        types[t++] = static_cast<int>(PyLong_AsLong(PyList_GetItem(pType,  i)));
        
        if(t==1000)
        {
            top->add_bonds(t, types, atom1, atom2);
            t = 0;
        }
    }
    
    if(t>0) top->add_bonds(t, types, atom1, atom2);
    return Py_BuildValue("i", top->nbonds);
}

PYAPI(add_angles)
{
    Topology *top;
    PyObject *pType, *pAtom1, *pAtom2, *pAtom3;
    PyArg_ParseTuple(args, "LOOOO", &top, &pType, &pAtom1, &pAtom2, &pAtom3);
    
    int t = 0, n = PyList_Size(pType);
    int types[1000], atom1[1000], atom2[1000], atom3[1000];
    
    for(int i=0; i<n; i++)
    {
        atom1[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom1, i)));
        atom2[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom2, i)));
        atom3[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom3, i)));
        types[t++] = static_cast<int>(PyLong_AsLong(PyList_GetItem(pType,  i)));
        
        if(t==1000)
        {
            top->add_angls(t, types, atom1, atom2, atom3);
            t = 0;
        }
    }
    
    if(t>0) top->add_angls(t, types, atom1, atom2, atom3);
    return Py_BuildValue("i", top->nangls);
}

PYAPI(add_dihedrals)
{
    Topology *top;
    PyObject *pType, *pAtom1, *pAtom2, *pAtom3, *pAtom4;
    PyArg_ParseTuple(args, "LOOOOO", &top, &pType, &pAtom1, &pAtom2, &pAtom3, &pAtom4);
    
    int t = 0, n = PyList_Size(pType);
    int types[1000], atom1[1000], atom2[1000], atom3[1000], atom4[1000];
    
    for(int i=0; i<n; i++)
    {
        atom1[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom1, i)));
        atom2[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom2, i)));
        atom3[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom3, i)));
        atom4[t]   = static_cast<int>(PyLong_AsLong(PyList_GetItem(pAtom4, i)));
        types[t++] = static_cast<int>(PyLong_AsLong(PyList_GetItem(pType,  i)));
        
        if(t==1000)
        {
            top->add_dihes(t, types, atom1, atom2, atom3, atom4);
            t = 0;
        }
    }
    
    if(t>0) top->add_dihes(t, types, atom1, atom2, atom3, atom4);
    return Py_BuildValue("i", top->ndihes);
}

PYAPI(build_special)
{
    Topology *top;
    int bond, angle, dihed;
    PyArg_ParseTuple(args, "Lppp", &top, &bond, &angle, &dihed);
    top->build_special(bond==1, angle==1, dihed==1);    
    Py_RETURN_NONE;
}
/*
PYAPI(get_pair_type)
{
    int i, j;
    PyArg_ParseTuple(args, "ii", &i, &j);
    return Py_BuildValue("i", get_pair_tid(i, j));
}
*/
static PyMethodDef cModPyMethods[] =
{
    {"create",     create,     METH_VARARGS, "Create topology."},
    {"destroy",    destroy,    METH_VARARGS, "Destroy topology."},
    {"add_atoms",  add_atoms,  METH_VARARGS, "Add atoms."},
    {"add_bonds",  add_bonds,  METH_VARARGS, "Add bonds."},
    {"add_angles", add_angles, METH_VARARGS, "Add angles."},
    {"add_dihedrals", add_dihedrals,  METH_VARARGS, "Add dihedrals."},
    {"build_special", build_special, METH_VARARGS, "Build special pair-lists."},
    //{"get_pair_type", get_pair_type, METH_VARARGS, "Return type id of a pairs."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_topol", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_topol(void)
{
    return PyModule_Create(&cModPy);
}
