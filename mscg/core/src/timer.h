#ifndef TIMER_H
#define TIMER_H

#include <ctime>

enum {IO, LIST, TABLE, MATRIX, SOLVER};

class Timer
{
  public:
    
    Timer();
    virtual ~Timer();
    
    int ntimers;
    int *count;
    double *seconds;
    
    clock_t start, last;
    
    void reset();
    void click(int);
    void report();
};

#endif