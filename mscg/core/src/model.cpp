#include "model.h"

Model::Model()
{
    this->tid = 0;
    this->list = 0;
    this->dF = 0;
    this->dU = 0;
    nparam  = 0;
}

Model::~Model()
{
    
}

void Model::setup(int tid, void* list, double *dF, double *dU)
{
    this->tid = tid;
    this->list = list;
    this->dF = dF;
    this->dU = dU;
}
