#include <Python.h>
#include <numpy/arrayobject.h>

#include "traj_trr.h"
#include "traj_lammps.h"
#include <string.h>

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETTRR() Traj *p; PyArg_ParseTuple(args, "L", &p)

template <typename T>
PyObject* open_traj(PyObject* self, PyObject* args)
{
    char *filename, *mode;
    PyArg_ParseTuple(args, "ss", &filename, &mode);
    T *f = new T(filename, mode);
    return Py_BuildValue("L", f);
}

PYAPI(close)
{
    GETTRR();
    delete p;
    Py_RETURN_NONE;
}

PyObject* rewind(PyObject* self, PyObject* args)
{
    GETTRR();
    p->rewind();
    Py_RETURN_NONE;
}

PyObject* get_status(PyObject* self, PyObject* args)
{
    GETTRR();
    PyObject *z = Py_BuildValue("i", p->status);
    return z;
}

PyObject* get_natoms(PyObject* self, PyObject* args)
{
    GETTRR();
    PyObject *z = Py_BuildValue("i", p->natoms);
    return z;
}

template<char ID>
PyObject* has_attr(PyObject* self, PyObject* args)
{
    GETTRR();
    int value = (ID=='t'?p->has_type:(ID=='v'?p->has_vel:p->has_force));
    PyObject *z = Py_BuildValue("O", value>0?Py_True:Py_False);
    return z;
}

PyObject* read_frame(PyObject* self, PyObject* args)
{
    Traj* traj;
    PyArrayObject *box, *type, *x, *v, *f;
    
    PyArg_ParseTuple(args, "LOOOOO", &traj, &box, &type, &x, &v, &f);
    int error = traj->read_next_frame();
    
    if(error == 0)
    {
        memcpy(box->data, traj->box, sizeof(float)*3);
        memcpy(x->data,   traj->x,   sizeof(float)*traj->natoms*3);
        
        if(traj->has_type)  memcpy(type->data, traj->type, sizeof(int)*traj->natoms);
        if(traj->has_vel)   memcpy(v->data,    traj->v,    sizeof(float)*traj->natoms*3);
        if(traj->has_force) memcpy(f->data,    traj->f,    sizeof(float)*traj->natoms*3);
    }
    
    PyObject *z = Py_BuildValue("O", error==0?Py_True:Py_False);
    return z;
}

PyObject* write_frame(PyObject* self, PyObject* args)
{
    Traj* traj;
    PyArrayObject *box, *type, *x, *v, *f;
    PyArg_ParseTuple(args, "LOOOOO", &traj, &box, &type, &x, &v, &f);
    
    traj->natoms = x->dimensions[0];
    traj->allocate();
    
    memcpy(traj->box, box->data, sizeof(float)*3);
    memcpy(traj->x,   x->data,   sizeof(float)*traj->natoms*3);
    
    traj->has_type  = Py_None != (PyObject *)type;
    traj->has_vel   = Py_None != (PyObject *)v;
    traj->has_force = Py_None != (PyObject *)f;
    
    if(traj->has_type)  memcpy(traj->type, type->data, sizeof(int)*traj->natoms);
    if(traj->has_vel)   memcpy(traj->v,    v->data,    sizeof(float)*traj->natoms*3);
    if(traj->has_force) memcpy(traj->f,    f->data,    sizeof(float)*traj->natoms*3);
    
    int error = traj->write_frame();
    PyObject *z = Py_BuildValue("O", error==0?Py_True:Py_False);
    return z;
}

static PyMethodDef cModPyMethods[] =
{
    {"open_trr",    open_traj<TrajTRR>,    METH_VARARGS, "Open Gromacs TRR file."},
    {"open_lmp",    open_traj<TrajLAMMPS>, METH_VARARGS, "Open LAMMPS cumstomized dump file."},
    {"close",       close,      METH_VARARGS, "Close the opened file."},
    {"rewind",      rewind,     METH_VARARGS, "Rewind the opened file."},
    {"get_status",  get_status, METH_VARARGS, "Get trajectory status."},
    {"get_natoms",  get_natoms, METH_VARARGS, "Return number of atoms."},
    {"has_type",    has_attr<'t'>, METH_VARARGS, "If trajectory has type data."},
    {"has_vel",     has_attr<'v'>, METH_VARARGS, "If trajectory has velocity data."},
    {"has_force",   has_attr<'f'>, METH_VARARGS, "If trajectory has force data."},
    {"read_frame",  read_frame, METH_VARARGS, "Read a frame."},
    {"write_frame", write_frame, METH_VARARGS, "Write a frame."},
    {NULL, NULL}
};

static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "cxx_traj", /* name of module */
    "",         /* module documentation, may be NULL */
    -1,         /* keeps state in global variables */
    cModPyMethods
};

PyMODINIT_FUNC PyInit_cxx_traj(void)
{
    import_array();
    return PyModule_Create(&cModPy);
}
