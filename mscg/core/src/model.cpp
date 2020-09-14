#include "model.h"

Model::Model(int tid, void* list, double *dF, double *dU)
{
    this->tid = tid;
    this->list = list;
    this->dF = dF;
    this->dU = dU;
    
    nparam  = 0;
}

Model::~Model()
{
    
}
