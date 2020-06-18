#include "py_model.h"
#include "model_pair_bspline.h"
#include "pair_list.h"

PYAPI(create)
{
    PairList *pair;
    int tid, order;
    double res, xmin; 
    PyArg_ParseTuple(args, "Liidd", &pair, &tid, &order, &res, &xmin);
    ModelPairBSpline *p = new ModelPairBSpline(pair, tid, order, res, xmin);
    return Py_BuildValue("L", p);
}

PYAPI(setup_cache)
{
    ModelPairBSpline *p;
    double ddx_factor;
    PyArg_ParseTuple(args, "Ld", &p, &ddx_factor);
    p->setup_cache(ddx_factor);
    Py_RETURN_NONE;
}
 
BEGIN_PY_API(ModelPairBSpline)
    DECLARE_API(create,      create,      "Create table object.")
    DECLARE_API(setup_cache, setup_cache, "Setup cache acceleration.")
END_PY_API()

DECLARE_PY_MODULE(cxx_model_pair_bspline)



    