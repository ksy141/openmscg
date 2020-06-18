#include "delta.h"

Delta::Delta()
{
    coeff = 0;
    ncols = 0;
}

Delta::~Delta()
{
    if(coeff) delete [] coeff;
}
