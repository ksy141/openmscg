#include "model.h"

Model::Model(void* list, int type_id)
{
    this->list = list;
    this->type_id = type_id;
    
    nparam  = 0;
    dudl    = 0;
}

Model::~Model()
{
    if(dudl) delete [] dudl;
}
