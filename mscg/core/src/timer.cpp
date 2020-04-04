#include "timer.h"
#include "std.h"

const char *timer_names[] = {
    "IO", "LIST", "TABLE", "MATRIX", "SOLVER",
    0
};

Timer::Timer()
{
    ntimers = 0;
    while(timer_names[ntimers]) ntimers++;
    
    count = new int[ntimers];
    seconds = new double[ntimers];
    
    for(int i=0; i<ntimers; i++)
    {
        count[i] = 0;
        seconds[i] = 0.0;
    }
}

Timer::~Timer()
{
    delete [] count;
    delete [] seconds;
}

void Timer::reset()
{
    last = clock();
    start = last;
}

void Timer::click(int id)
{
    clock_t now = clock();
    
    count[id]++;
    seconds[id] += (double)(now - last)/CLOCKS_PER_SEC;
    
    last = now;
}

void Timer::report()
{
    clock_t now = clock();
    
    printf("Total time: %0.2lf seconds\n", (double)(now - start)/CLOCKS_PER_SEC);
    
    for(int i=0; i<ntimers; i++)
        printf("%8s %6.2lf (%d)\n", timer_names[i], seconds[i], count[i]);
}