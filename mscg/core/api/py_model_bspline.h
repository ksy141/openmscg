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
    int order;

    PyArg_ParseTuple(args, "dddi", &xmin, &xmax, &resolution, &order);
    MODEL_CLASS *p = new MODEL_CLASS(xmin, xmax, resolution, order);
    return Py_BuildValue("L", p);
}

PYAPI(setup)
{
    MODEL_CLASS *p;
    int tid;
    void *plist;
    PyArrayObject *dF, *dU;

    PyArg_ParseTuple(args, "LiLOO", &p, &tid, &plist, &dF, &dU);
    p->setup(tid, plist, NP_DATA(dF), NP_DATA(dU));
    Py_RETURN_NONE;
}

PYAPI(setup_cache)
{
    MODEL_CLASS *p;
    double ddx_factor;
    PyArg_ParseTuple(args, "Ld", &p, &ddx_factor);
    p->setup_cache(ddx_factor);
    Py_RETURN_NONE;
}

BEGIN_PY_API(MODEL_CLASS)
    DECLARE_API(get_npars,   get_npars,   "Get count of parameters.")
    DECLARE_API(create,      create,      "Create table object.")
    DECLARE_API(setup,       setup,       "Setup model attributes.")
    DECLARE_API(setup_cache, setup_cache, "Setup cache acceleration.")
#ifdef PY_API_EXTRA
    PY_API_EXTRA
#endif
END_PY_API()
