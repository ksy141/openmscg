#include "py_model.h"
#include "model_pair_bspline.h"
#include "pair_list.h"

PYAPI(get_npars)
{
    double xmin, xmax, resolution;
    int order;
    
    PyArg_ParseTuple(args, "dddi", &xmin, &xmax, &resolution, &order);
    return Py_BuildValue("i", BSpline::get_nbreak(xmin, xmax, resolution) + order - 2);
}

PYAPI(create)
{
    double xmin, xmax, resolution;
    int order, tid;
    void *plist;
    PyArrayObject *dF, *dU;
    
    PyArg_ParseTuple(args, "dddiiLOO", &xmin, &xmax, &resolution, &order, &tid, &plist, &dF, &dU);
    ModelPairBSpline *p = new ModelPairBSpline(xmin, xmax, resolution, order, tid, plist, NP_DATA(dF), NP_DATA(dU));
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
    DECLARE_API(get_npars,   get_npars,   "Get count of parameters.")
    DECLARE_API(create,      create,      "Create table object.")
    DECLARE_API(setup_cache, setup_cache, "Setup cache acceleration.")
END_PY_API()

DECLARE_PY_MODULE(cxx_model_pair_bspline)

    
    