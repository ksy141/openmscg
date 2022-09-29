#include "py_model.h"
#include "model_nb3b_bspline.h"


#define MODEL_CLASS ModelNB3BBSpline

PYAPI(setup_ex)
{
    MODEL_CLASS *p;
    int tid_ij, tid_ik;
    double gamma_ij, gamma_ik;
    double a_ij, a_ik;

    PyArg_ParseTuple(args, "Liddidd", &p, &tid_ij, &gamma_ij, &a_ij,
        &tid_ik, &gamma_ik, &a_ik);

    p->setup_ex(tid_ij, gamma_ij, a_ij, tid_ik, gamma_ik, a_ik);
    Py_RETURN_NONE;
}

#define PY_API_EXTRA {"setup_ex", setup_ex, METH_VARARGS, "Extra parameters."},

#include "py_model_bspline.h"
DECLARE_PY_MODULE(cxx_model_nb3b_bspline)

#undef PY_API_EXTRA
