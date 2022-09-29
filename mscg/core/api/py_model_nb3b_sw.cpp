#include "py_model.h"
#include "model_nb3b_sw.h"

PYAPI(create)
{
    double gamma_ij, gamma_ik, a_ij, a_ik, theta0;

    PyArg_ParseTuple(args, "ddddd", &gamma_ij, &a_ij, &gamma_ik, &a_ik, &theta0);
    ModelNB3BSW *p = new ModelNB3BSW(gamma_ij, a_ij, gamma_ik, a_ik, theta0);
    return Py_BuildValue("L", p);   
}

PYAPI(setup)
{
    ModelNB3BSW *p;
    int tid;
    void *plist;
    PyArrayObject *dF, *dU;

    PyArg_ParseTuple(args, "LiLOO", &p, &tid, &plist, &dF, &dU);
    p->setup(tid, plist, NP_DATA(dF), NP_DATA(dU));
    Py_RETURN_NONE;
}


PYAPI(setup_ex)
{
    ModelNB3BSW *p;
    int tid_ij, tid_ik;

    PyArg_ParseTuple(args, "Lii", &p, &tid_ij, &tid_ik);
    p->setup_ex(tid_ij, tid_ik);
    Py_RETURN_NONE;
}


BEGIN_PY_API(ModelNB3BSW)
  DECLARE_API(create, create, "Create model object.")
  DECLARE_API(setup, setup, "Setup parameters.")
  DECLARE_API(setup_ex, setup_ex, "Setup extra parameters.")
END_PY_API()

DECLARE_PY_MODULE(cxx_model_nb3b_sw)
