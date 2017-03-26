#include <iostream>
#include <iomanip>
#include <string.h>
#include <fstream>
#include <cmath>
#include <vector>

#define pi 3.14159265358979323846 

using namespace std;

extern "C" void dstev_(char* job, int* N, double* D, double* OFFD, double* EV, int* VDIM, double* WORK, int* INFO);
extern "C" void dsyev_(char* JOB, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO);

void print_to_file(const char*, vector<double> &, vector<double> &, int);
void square_guide_setup(vector<double> &, int, double, double, double);
void fd_matrix(vector<double> &, vector<double> &, int, double, double, double);

int main(){
    
    /* PROGRAM PARAMETERS */ 
    int dim = 40, i = 0, j=0;
    int n = dim*dim;
    double k = pi/1.55e-6;
    double step = 1e-7;
    vector<double> eps(n, 0);
    char vectorfile[30];

    /* LAPACK PARAMETERS */
    char job = 'V';
    char uplo = 'U';
    int N = n;
    vector<double> a(n*n, 0);
    int lda = n;
    vector<double> w(n, 0);
    int lwork = 3*n -1;
    vector<double> work(lwork, 0);
    int info;


    
    double ratio = 1/25.0;
    double n_out = 1.4, n_in = 1.5;
    
    square_guide_setup(eps, dim, ratio, n_out, n_in);
    fd_matrix(a, eps, dim, step, step, k);

    dsyev_(&job, &uplo, &n, &a[0], &lda, &w[0], &work[0], &lwork, &info);
    if(info == 0){
        cout << "Solution Found" << endl;
    }else{
        cout << "LAPACK Failed\nFORTRAN routine exited with INFO = "<< info << endl;
    }    

/*
    for(i=0; i<dim*dim; i++){
        if(i%dim == 0){cout << endl;}       
        cout << eps[i] << "\t";
    }

    cout << "\n\n\n\n\n";
    cout << setprecision(2);
    for(i=0; i<n*n; i++){
        if(i%n == 0){cout << endl;}       
        cout << a[i] << "\t";
    }*/

    cout << "Printing first Eigenvector to file" << endl;
    vector<double> efield(a.begin() + n*(n-1), a.begin() + n*(n));
    vector<double> distance;
    for(i=int(-n/2); i<int(n/2); i++){
        distance.push_back(i*step);
    }

    sprintf(vectorfile, "./TM%d.txt", 0);
    print_to_file(vectorfile, distance, efield, n);

    return 0;
}

void print_to_file(const char* filename, vector<double> &var1, vector<double> &var2, int size){
    /* 
     * Prints data stored in var1 and var2 in two column format in the 
     * file called 'filename', where the length of the arrays var1/2 must
     * be the same.
     * Size is the size of the arrays var1/2 
     * */
    ofstream myfile(filename);
    int i = 0;
    if(var2.empty()){
        int total = var1.size();
        for(i=total-1; i>=total - size; i--){
            myfile << var1[i] << endl;
        }
    }
    else{
        for(i=0; i<size; i++){
            myfile << var1[i] << "\t" << var2[i] << endl;
        }
    }
    myfile.close();

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


void fd_matrix(vector<double> &a, vector<double> &eps, int dim, double dx, double dy, double k){
/*
 * IN:
 *      a -> vector of dimension dim^4 filled with just zeros (represents FD matrix)
 *      eps -> matrix containing epsilon_r at each coordinate (x, y)
 *      dim -> dimension in of of 1D problem
 *      dx -> step in x direction
 *      dy -> step in y direction
 *      k -> wave vector
 *
 * OUT:
 *      a -> vector now filled with elements of FD matrix to be passed into the lapack routine
 *           stored as a vector in order to be compatible with the LAPACK routine
 */
    int n = dim * dim;
    int i = 0, j=0;
    double a1 = 1/(dx*dx);
    double a2 = 1/(dx*dx);
    vector<double> a3(n*n, 0);
    vector<double> a4(n*n, 0);
    vector<double> a5(n*n, 0);

    /* INITIALIZING MULTIPLIERS */
    for(i=0; i<n*n; i++){
        if(i == 0){
            a3[i] = 1/(dy*dy);
            a4[i] = 1/(dy*dy) * 2 * eps[(i+1)%n]/(eps[i%n] + eps[(i+1)%n]);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k * eps[i%n];
        }
        else if(i ==n*n-1){
            a3[i] = 1/(dy*dy) * 2 * eps[(i-1)%n]/(eps[i%n] + eps[(i-1)%n]);
            a4[i] = 1/(dy*dy);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k * eps[i%n];
        }
        else{
            a3[i] = 1/(dy*dy) * 2 * eps[(i-1)%n]/(eps[i%n] + eps[(i-1)%n]);
            a4[i] = 1/(dy*dy) * 2 * eps[(i+1)%n]/(eps[i%n] + eps[(i+1)%n]);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k * eps[i%n];
        }
    }

    
    /* CONTRUCTING MATRIX */
    for(i=0; i < n; i++){
        if( i == 0){
            a[0] = a5[i];
            a[1] = a4[i*n + 1];
            a[dim] = a2;
        }
        else if(i == n -1){
            a[i*n + i] = a5[i*n + i];
            a[i*n + i-1] = a3[i*n + i -1];
            a[i*n + i - dim] = a3[i*n + i - dim];
        }else if(i < dim){
            a[i*n + i -1] = a3[i*n + i -1];
            a[i*n + i] = a5[i*n + i];
            a[i*n + i+1] = a4[i*n + i+1];
            a[i*n + i + dim] = a2;            
        }else if(i >= n - dim){
            a[i*n + i - dim] = a1;
            a[i*n + i -1] = a3[i*n + i -1];
            a[i*n + i] = a5[i*n + i];
            a[i*n + i+1] = a4[i*n + i+1];
        }else{
            a[i*n + i - dim] = a1;
            a[i*n + i -1] = a3[i*n + i -1];
            a[i*n + i] = a5[i*n + i];
            a[i*n + i+1] = a4[i*n + i+1];
            a[i*n + i + dim] = a2; 
        }
    }
}  
