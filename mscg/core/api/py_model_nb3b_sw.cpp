#include "py_model.h"
#include "model_nb3b_sw.h"


#define MODEL_CLASS ModelNB3BSW

PYAPI(create)
{
    int tid_ij, tid_ik;
    double gamma_ij, gamma_ik, a_ij, a_ik, theta0;

    PyArg_ParseTuple(args, "iddiddd", &tid_ij, &gamma_ij, &a_ij,
       &tid_ik, &gamma_ik, &a_ik, &theta0);
    
    
    MODEL_CLASS *p = new ModelNB3BSW(tid_ij, gamma_ij, a_ij, tid_ik, gamma_ik, a_ik, theta0);
    return Py_BuildValue("L", p);   
}

BEGIN_PY_API(MODEL_CLASS)
  DECLARE_API("create", create, "Create model object.")
END_PY_API()

DECLARE_PY_MODULE(cxx_model_nb3b_sw)
