#include "matrix.h"
#include "traj_trr.h"
#include "pair_list.h"
#include "table_pair_bspline.h"
#include "timer.h"

#include "std.h"

int main()
{
    const char *filename = "CGTraj.trr";
    
    Estimator *estimator = new Estimator();
    
    TrajTRR *trr = new TrajTRR(filename);
    
    estimator->set_natoms(trr->natoms);
    
    PairList *pair = new PairList(trr->natoms, 10.0, 5.0);
    
    TablePairBSpline * pair_bs = new TablePairBSpline(estimator, pair, 3, 0.1, 2.7);
    pair_bs->setup_table();
    
    estimator->setup();
    
    Timer *tm = new Timer();
    tm->reset();
    
    for(int i=0; i<1000; i++)
    {
        trr->read_next_frame();
        
        for(int i=0; i<trr->natoms; i++) for(int d=0; d<3; d++)
            if(trr->x[i][d]<0) trr->x[i][d] += trr->box[d];
            else if(trr->x[i][d]>trr->box[d]) trr->x[i][d] -= trr->box[d];
        
        tm->click(IO);
        
        pair->setup_bins(trr->box);
        pair->build(trr->x);
        
        tm->click(LIST);
        
        estimator->empty_coeff_matrix();
        
        pair_bs->process();
        
        tm->click(TABLE);
        
        estimator->process_coeff_matrix((float*)trr->f);
        
        tm->click(MATRIX);
    }
    
    estimator->solve();
    
    tm->click(SOLVER);
    
    tm->report();
    
    print_vector(estimator->vector_cov, estimator->ncols);

    double *output = new double[100];
    pair_bs->dump(estimator->vector_cov, output, 0.05, 100);
    print_vector(output, 100);  
}
