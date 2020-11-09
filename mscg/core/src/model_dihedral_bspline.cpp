#include "model_dihedral_bspline.h"
#include "bond_list.h"
#include "defs.h"

ModelDihedralBSpline::ModelDihedralBSpline(double xmin, double xmax, double resolution, int order) : BSpline(order, resolution, xmin, xmax)
{
    nparam = ncoeff;
}

ModelDihedralBSpline::~ModelDihedralBSpline()
{

}

void ModelDihedralBSpline::compute_fm()
{
    BondList *lst = (BondList*)(this->list);
    vec4i *atoms = lst->dihedral_atoms;
    int* types = lst->dihedral_types;
    float *phi = lst->phi_dihedral;
    
    float *dpd1x = lst->dpd1x_dihedral;
    float *dpd1y = lst->dpd1y_dihedral;
    float *dpd1z = lst->dpd1z_dihedral;
    float *dpd2x = lst->dpd2x_dihedral;
    float *dpd2y = lst->dpd2y_dihedral;
    float *dpd2z = lst->dpd2z_dihedral;
    float *dpd3x = lst->dpd3x_dihedral;
    float *dpd3y = lst->dpd3y_dihedral;
    float *dpd3z = lst->dpd3z_dihedral;
    float *dpd4x = lst->dpd4x_dihedral;
    float *dpd4y = lst->dpd4y_dihedral;
    float *dpd4z = lst->dpd4z_dihedral;
    
    for(int i=0; i<lst->ndihedrals; i++) if(types[i] == tid)
    {         
        double *b;
        size_t istart;
        int nn;
        eval_coeffs(phi[i], &b, &istart, &nn);

        int i1 = atoms[i][0], i2 = atoms[i][1], i3 = atoms[i][2], i4 = atoms[i][3];
        double *coeff_1[3], *coeff_2[3], *coeff_3[3], *coeff_4[3];

        coeff_1[0] = dF + i1 * 3 * nparam;
        coeff_2[0] = dF + i2 * 3 * nparam;
        coeff_3[0] = dF + i3 * 3 * nparam;
        coeff_4[0] = dF + i4 * 3 * nparam;
        
        coeff_1[1] = dF + (i1 * 3 + 1) * nparam;
        coeff_2[1] = dF + (i2 * 3 + 1) * nparam;
        coeff_3[1] = dF + (i3 * 3 + 1) * nparam;
        coeff_4[1] = dF + (i4 * 3 + 1) * nparam;
        
        coeff_1[2] = dF + (i1 * 3 + 2) * nparam;
        coeff_2[2] = dF + (i2 * 3 + 2) * nparam;
        coeff_3[2] = dF + (i3 * 3 + 2) * nparam;
        coeff_4[2] = dF + (i4 * 3 + 2) * nparam;
        
        for(int c=0; c<nn; c++)
        {
            double Bi = b[c];
            int pos = istart + c;
            
            coeff_1[0][pos] += Bi * dpd1x[i];
            coeff_1[1][pos] += Bi * dpd1y[i];
            coeff_1[2][pos] += Bi * dpd1z[i];
            
            coeff_2[0][pos] += Bi * dpd2x[i];
            coeff_2[1][pos] += Bi * dpd2y[i];
            coeff_2[2][pos] += Bi * dpd2z[i];
            
            coeff_3[0][pos] += Bi * dpd3x[i];
            coeff_3[1][pos] += Bi * dpd3y[i];
            coeff_3[2][pos] += Bi * dpd3z[i];
            
            coeff_4[0][pos] += Bi * dpd4x[i];
            coeff_4[1][pos] += Bi * dpd4y[i];
            coeff_4[2][pos] += Bi * dpd4z[i];
        }
    }
}

void ModelDihedralBSpline::compute_rem()
{
    
}

void ModelDihedralBSpline::get_table(double *params, double* in, double* out, int size)
{
    eval(params, in, out, size);
}
