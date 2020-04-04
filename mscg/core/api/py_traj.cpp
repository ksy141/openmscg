#include <Python.h>
#include "traj_trr.h"
#include "traj_lammps.h"

#define PYAPI(api) PyObject* api(PyObject* self, PyObject* args)
#define GETTRR() Traj *p; PyArg_ParseTuple(args, "L", &p)

PyObject* open_trr(PyObject* self, PyObject* args)
{
    char *filename;
    PyArg_ParseTuple(args, "s", &filename);
    TrajTRR *f = new TrajTRR(filename);
    return Py_BuildValue("L", f);
}

PyObject* open_lmp(PyObject* self, PyObject* args)
{
    char *filename;
    PyArg_ParseTuple(args, "s", &filename);
    TrajLAMMPS *f = new TrajLAMMPS(filename);
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

PyObject* has_force(PyObject* self, PyObject* args)
{
    GETTRR();
    PyObject *z = Py_BuildValue("O", p->has_force>0?Py_True:Py_False);
    return z;
}

PyObject* get_box(PyObject* self, PyObject* args)
{
    GETTRR();
    PyObject *z = Py_BuildValue("[f,f,f]", p->box[0], p->box[1], p->box[2]);
    return z;
}

PyObject* next_frame(PyObject* self, PyObject* args)
{
    GETTRR();
    int error = p->read_next_frame();
    PyObject *z = Py_BuildValue("O", error==0?Py_True:Py_False);
    return z;
}

PyObject* get_vector(PyObject* args, int which)
{
    GETTRR();
    int nn = p->natoms * 3;
    float *src = (which==0?((float*)(p->x)):((float*)(p->f)));
    PyObject *des  = PyList_New(nn);
    
    for(int i=0; i<nn; i++) PyList_SetItem(des, i, Py_BuildValue("f", src[i]));
    return des;
}

PyObject* get_x(PyObject* self, PyObject* args) { return get_vector(args, 0); }
PyObject* get_f(PyObject* self, PyObject* args) { return get_vector(args, 1); }

static PyMethodDef cModPyMethods[] =
{
    {"open_trr",   open_trr,   METH_VARARGS, "Open Gromacs TRR file."},
    {"open_lmp",   open_lmp,   METH_VARARGS, "Open LAMMPS cumstomized dump file."},
    {"close",      close,      METH_VARARGS, "Close the opened file."},
    {"rewind",     rewind,     METH_VARARGS, "Rewind the opened file."},
    {"get_status", get_status, METH_VARARGS, "Get trajectory status."},
    {"get_natoms", get_natoms, METH_VARARGS, "Return number of atoms."},
    {"has_force",  has_force,  METH_VARARGS, "If trajectory has force data."},
    {"get_box",    get_box,    METH_VARARGS, "Return box size."},
    {"get_x",      get_x,      METH_VARARGS, "Return x."},
    {"get_f",      get_f,      METH_VARARGS, "Return f."},
    {"next_frame", next_frame, METH_VARARGS, "Read a frame."},
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
    return PyModule_Create(&cModPy);
}
