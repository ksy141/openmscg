#include <Python.h>
#include <numpy/arrayobject.h>

#define PYAPI(api) static PyObject* api(PyObject* self, PyObject* args)

template<class T>
class PyModel
{
public:
    
    PYAPI(Test)
    {
        int value;
        PyArg_ParseTuple(args, "i", &value);
        printf("Testing Value = %d\n", value);
        Py_RETURN_NONE;
    }
    
    PYAPI(Destroy)
    {
        T *p;
        PyArg_ParseTuple(args, "L", &p);
        delete p;
        Py_RETURN_NONE;
    }
        
    PYAPI(ComputeDFDL)
    {
        T *p;
        PyArg_ParseTuple(args, "L", &p);
        p->compute_dfdl();
        Py_RETURN_NONE;
    }
    
    
    PYAPI(ComputeDUDL)
    {
        T *p;
        PyArg_ParseTuple(args, "L", &p);
        p->compute_dudl();

        PyObject *dudl = PyList_New(p->nparam);
        for(int i=0; i<p->nparam; i++) PyList_SetItem(dudl, i, Py_BuildValue("f", p->dudl[i]));
        return dudl;
    }
    
    PYAPI(ComputeEnergyTable)
    {
        T *p;
        PyArrayObject *pars, *in, *out;
        
        PyArg_ParseTuple(args, "LOOO", &p, &pars, &in, &out);
        p->compute_etable((double*)(pars->data), (double*)(in->data), (double*)(out->data), in->dimensions[0]);
        Py_RETURN_NONE;
    }
};

#define DECLARE_API(name, method, desc) {#name, method, METH_VARARGS, desc},

#define BEGIN_PY_API(class) \
static PyMethodDef cModPyMethods[] = {\
    {"test", PyModel<class>::Test, METH_VARARGS, "Just for test."}, \
    {"destroy", PyModel<class>::Destroy, METH_VARARGS, "Destroy model."}, \
    {"compute_dfdl", PyModel<class>::ComputeDFDL, METH_VARARGS, "Compute dF/dLambda."}, \
    {"compute_dudl", PyModel<class>::ComputeDUDL, METH_VARARGS, "Compute dU/dLambda."}, \
    {"compute_etable", PyModel<class>::ComputeEnergyTable, METH_VARARGS, "Compute energy table."},

#define END_PY_API() {NULL, NULL}};

#define DECLARE_PY_MODULE(name) static struct PyModuleDef cModPy = { \
    PyModuleDef_HEAD_INIT, #name, "", -1, cModPyMethods }; \
    PyMODINIT_FUNC PyInit_##name(void) { \
        import_array(); \
        return PyModule_Create(&cModPy); }
