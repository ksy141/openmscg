#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "traj_lammps.h"
#include "traj_trr.h"
#include "traj_dcd.h"
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

PyObject* get_timestep(PyObject* self, PyObject* args)
{
    GETTRR();
    PyObject *z = Py_BuildValue("i", p->step);
    return z;
}

PyObject* set_timestep(PyObject* self, PyObject* args)
{
    Traj* traj;
    int ts;
    PyArg_ParseTuple(args, "Li", &traj, &ts);
    traj->step = ts;
    Py_RETURN_NONE;
}

PyObject* has_attr(PyObject* self, PyObject* args)
{
    Traj* traj;
    char* attr;
    PyArg_ParseTuple(args, "Ls", &traj, &attr);
    return Py_BuildValue("O", traj->attrs[attr[0]]>0?Py_True:Py_False);
}

PyObject* read_frame(PyObject* self, PyObject* args)
{
    GETTRR();
    int error = p->read_next_frame();
    return Py_BuildValue("O", error==0?Py_True:Py_False);
}

#define NOT_NONE(p) (Py_None != (PyObject *)p)

PyObject* get_frame(PyObject* self, PyObject* args)
{    
    Traj* traj;
    PyArrayObject *box, *t, *q, *x, *v, *f;
    PyArg_ParseTuple(args, "LOOOOOO", &traj, &box, &t, &q, &x, &v, &f);
    
    memcpy(PyArray_DATA(box), traj->box, sizeof(float)*3);
    memcpy(PyArray_DATA(x),   traj->x,   sizeof(float)*traj->natoms*3);
    
    if NOT_NONE(t) memcpy(PyArray_DATA(t), traj->t, sizeof(int)*traj->natoms);
    if NOT_NONE(q) memcpy(PyArray_DATA(q), traj->q, sizeof(float)*traj->natoms);
    if NOT_NONE(v) memcpy(PyArray_DATA(v), traj->v, sizeof(float)*traj->natoms*3);
    if NOT_NONE(f) memcpy(PyArray_DATA(f), traj->f, sizeof(float)*traj->natoms*3);
    
    Py_RETURN_NONE;
}

PyObject* write_frame(PyObject* self, PyObject* args)
{
    Traj* traj;
    PyArrayObject *box, *t, *q, *x, *v, *f;
    PyArg_ParseTuple(args, "LOOOOOO", &traj, &box, &t, &q, &x, &v, &f);
    
    traj->natoms = PyArray_DIMS(x)[0];
    traj->allocate();
    
    memcpy(traj->box, PyArray_DATA(box), sizeof(float)*3);
    memcpy(traj->x,   PyArray_DATA(x),   sizeof(float)*traj->natoms*3);
    
    traj->attrs['t'] = NOT_NONE(t);
    traj->attrs['q'] = NOT_NONE(q);
    traj->attrs['v'] = NOT_NONE(v);
    traj->attrs['f'] = NOT_NONE(f);
    
    if(traj->attrs['t']) memcpy(traj->t, PyArray_DATA(t), sizeof(int)*traj->natoms);
    if(traj->attrs['q']) memcpy(traj->q, PyArray_DATA(q), sizeof(float)*traj->natoms);
    if(traj->attrs['v']) memcpy(traj->v, PyArray_DATA(v), sizeof(float)*traj->natoms*3);
    if(traj->attrs['f']) memcpy(traj->f, PyArray_DATA(f), sizeof(float)*traj->natoms*3);
    
    int error = traj->write_frame();
    PyObject *z = Py_BuildValue("O", error==0?Py_True:Py_False);
    return z;
}

static PyMethodDef cModPyMethods[] =
{
    {"open_lmp",    open_traj<TrajLAMMPS>, METH_VARARGS, "Open LAMMPS cumstomized dump file."},
    {"open_trr",    open_traj<TrajTRR>,    METH_VARARGS, "Open Gromacs TRR file."},
    {"open_dcd",    open_traj<TrajDCD>,    METH_VARARGS, "Open CHARMM/NAMD DCD file."},
    {"close",       close,       METH_VARARGS, "Close the opened file."},
    {"rewind",      rewind,      METH_VARARGS, "Rewind the opened file."},
    {"get_status",  get_status,  METH_VARARGS, "Get trajectory status."},
    {"get_natoms",  get_natoms,  METH_VARARGS, "Return number of atoms."},
    {"get_timestep",get_timestep,METH_VARARGS, "Return timestep."},
    {"set_timestep",set_timestep,METH_VARARGS, "Set timestep."},
    {"has_attr",    has_attr,    METH_VARARGS, "If trajectory has the attribute."},
    {"read_frame",  read_frame,  METH_VARARGS, "Read a frame."},
    {"get_frame",   get_frame,   METH_VARARGS, "Get a copy of frame data."},
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
