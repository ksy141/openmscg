#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#define PYAPI(api) static PyObject* api(PyObject* self, PyObject* args)
#define NP_DATA(p) ((Py_None==(PyObject *)p)?0:(double*)PyArray_DATA(p))

template<class T>
class PyModel
{
public:
    
    PYAPI(Destroy)
    {
        T *p;
        PyArg_ParseTuple(args, "L", &p);
        delete p;
        Py_RETURN_NONE;
    }
    
    PYAPI(ComputeFM)
    {
        T *p;
        PyArg_ParseTuple(args, "L", &p);
        p->compute_fm();
        Py_RETURN_NONE;
    }
    
    PYAPI(ComputeREM)
    {
        T *p;
        PyArg_ParseTuple(args, "L", &p);
        p->compute_rem();
        Py_RETURN_NONE;
    }
    
    PYAPI(GetTable)
    {
        T *p;
        PyArrayObject *pars, *in, *out;
        
        PyArg_ParseTuple(args, "LOOO", &p, &pars, &in, &out);
        //printf("ABC %d %d\n", PyArray_DIMS(pars)[0], p->ncoeff);
        p->get_table((double*)PyArray_DATA(pars), (double*)PyArray_DATA(in), (double*)PyArray_DATA(out), PyArray_DIMS(in)[0]);
        Py_RETURN_NONE;
    }
};

#define DECLARE_API(name, method, desc) {#name, method, METH_VARARGS, desc},

#define BEGIN_PY_API(class) \
static PyMethodDef cModPyMethods[] = {\
    {"destroy", PyModel<class>::Destroy, METH_VARARGS, "Destroy model."}, \
    {"compute_rem", PyModel<class>::ComputeREM, METH_VARARGS, "Compute dF/dLambda for REM."}, \
    {"compute_fm", PyModel<class>::ComputeFM, METH_VARARGS, "Compute dU/dLambda for FM."}, \
    {"get_table", PyModel<class>::GetTable, METH_VARARGS, "Return the model table."},

#define END_PY_API() {NULL, NULL}};

#define DECLARE_PY_MODULE(name) static struct PyModuleDef cModPy = { \
    PyModuleDef_HEAD_INIT, #name, "", -1, cModPyMethods }; \
    PyMODINIT_FUNC PyInit_##name(void) { \
        import_array(); \
        return PyModule_Create(&cModPy); }