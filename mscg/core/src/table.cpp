#include "table.h"

Table::Table()
{
    coeff = 0;
    ncols = 0;
}

Table::~Table()
{
    if(coeff) delete [] coeff;
}
