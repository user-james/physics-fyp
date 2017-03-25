#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <vector>

#define pi 3.14159265358979323846 

using namespace std;

extern "C" void dstev_(char* job, int* N, double* D, double* OFFD, double* EV, int* VDIM, double* WORK, int* INFO);

void square_guide_setup(vector<double> &, int, double, double, double);


int main(){
    

    int dim = 4, i = 0, j=0;
    
     
    cout << "Enter size of square matrix: ";
    cin >> dim;
    vector<double> eps(dim*dim + 2, 0);

    double ratio = 1/3.0;
    double n_out = 1.4, n_in = 1.5;

    square_guide_setup(eps, dim, ratio, n_out, n_in);

    for(i=0; i<dim*dim; i++){
        if(i%dim == 0){cout << endl;}       
        cout << eps[i] << "\t";
    }

    cout << endl;

    return 0;
}



void square_guide_setup(vector<double> &eps , int dim, double centre_ratio, double n_out, double n_in){
/*
 * This function takes in an empty matrix and fills it with a symmetric square wave guide based on the size
 * of ratio you want in the wave guide
 *
 * Note: always makes a symmetric square waveguide but it may not be the correct ratio is the dimension 
 * passed in doesn't allow for a symmetric square wave guide of that ratio
 *
 */ 
    int left, right, i = 0, j = 0;
    left = dim * (1. - centre_ratio) / 2;
    right = dim  - 1 - left;

    for(i=0; i<dim; i++){
        for(j = 0; j<dim; j++){
            if( i<left || j<left){
                eps[i*dim + j] = n_out;
            }else if(i>right || j>right){
                eps[i*dim + j] = n_out;
            }else{
                eps[i*dim + j] = n_in;
            }            
        }
    }
 
}


void fd_matrix(double *a, vector< vector<double> > &eps, int dim, double dx, double dy, double k){
/*
 * IN:
 *      a -> vector of dimension dim^4 filled with just zeros (represents FD matrix)
 *      eps -> matrix containing epsilon_r at each coordinate (x, y)
 *      dim -> dimension in of of 1D problem
 *      dx -> step in x direction
 *      dy -> step in y direction
 *      k -> wave vector
 *
 * MODIFIED:
 *      a -> vector now filled with elements of FD matrix to be passed into the lapack routine
 *           stored as a vector in order to be compatible with the LAPACK routine
 */
    int i = 0, j=0;
    double a1 = 1/(dx*dx);
    double a2 = 1/(dx*dx);
    double a3[dim*dim];
    double a4[dim*dim];
    double a5[dim*dim];


    /* INITIALIZING MULTIPLIERS */
    for(i=0; i<dim*dim; i++){
        if(i == 0){
            a3[i] = 1/(dy*dy);
            a4[i] = 1/(dy*dy) * 2 * eps[i+1]/(eps[i] + eps[i+1]);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k * eps[i];
        }
        else if(i ==dim*dim-1){
            a3[i] = 1/(dy*dy) * 2 * eps[i-1]/(eps[i] + eps[i-1]);
            a4[i] = 1/(dy*dy);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k * eps[i];
        }
        else{
            a3[i] = 1/(dy*dy) * 2 * eps[i-1]/(eps[i] + eps[i-1]);
            a4[i] = 1/(dy*dy) * 2 * eps[i+1]/(eps[i] + eps[i+1]);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k * eps[i];
        }
    }

    /* CONTRUCTING MATRIX */
    for(i=0; i < dim*dim; i++){
        for(j=0; j<dim; j++){
            if(i >= dim && i < dim*(dim -1)){
            
            }
        }
    }
}  
