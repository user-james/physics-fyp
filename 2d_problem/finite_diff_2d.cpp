#include <iostream>
#include <iomanip>
#include <string.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

#define pi 3.14159265358979323846 

using namespace std;

extern "C" void dstev_(char* job, int* N, double* D, double* OFFD, double* EV, int* VDIM, double* WORK, int* INFO);
extern "C" void dsyev_(char* JOB, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO);
extern "C" void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA, double* WR, double* WI, double* VL, int* LDVL, double* VR, int* LDVR, double* WORK, int* LWORK, int* INFO);

void print_to_file_3d(const char*, vector<double> &, vector<double> &, vector<double> &, int);
void square_guide_setup(vector<double> &, int, double, double, double);
void fd_matrix(vector<double> &, vector<double> &, int, double, double, double);

int main(){
    
    /* PROGRAM PARAMETERS */ 
    int dim = 100, i = 0, j=0;
    int n = dim*dim;
    double k = pi/1.55e-6;
    double step = 8e-6/5;
    vector<double> eps(n, 0);
    char vectorfile[50];

    /* LAPACK PARAMETERS */
    char jobvl = 'N';
    char jobvr = 'V';
    int N = n;
    int lda = n;
    int ldvl = 1;
    int ldvr = n;
    int lwork = 6*n;          // chosen optimally after test runs
    int info;
    vector<double> wr(n);
    vector<double> wi(n);
    vector<double> a(lda*n, 0);
    vector<double> vl(ldvl*n, 0);
    vector<double> vr(ldvr*n, 0);
    vector<double> work(lwork, 0);


    
    double ratio = 1/9.0;
    double n_out = 3.16, n_in = 3.5;

    cout << "Constructing FD Matrix" << endl;    
    square_guide_setup(eps, dim, ratio, n_out, n_in);
    fd_matrix(a, eps, dim, step, step, k);

/*    
    for(i=0; i<dim*dim; i++){
        if(i%dim == 0){cout << endl;}       
        cout << eps[i] << "\t";
    }

    cout << "\n\n\n\n\n";
    cout << setprecision(3);
    for(i=0; i<n*n; i++){
        if(i%n == 0){cout << endl;}       
        cout << a[i] << "\t";
    }
    cout << endl;
*/
    
    cout << "Solving System" << endl;
    const clock_t begin_time = clock();
    
    /* SOLVING SYSTEM */
    dgeev_(&jobvl, &jobvr, &n, &a[0], &lda, &wr[0], &wi[0], &vl[0], &ldvl, &vr[0], &ldvr, &work[0], &lwork, &info);

    //dsyev_(&job, &uplo, &n, &a[0], &lda, &w[0], &work[0], &lwork, &info);
    if(info == 0){
        cout << "\nSolution Found" << endl;
        cout << "Time Taken: " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
        cout << "Optimal LWORK  = " << work[0] << endl;
    }else{
        cout << "LAPACK Failed\nFORTRAN routine exited with INFO = "<< info << endl;
        cout << "Time Taken: " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
    }    



    /* CREATING POINTS ON X-Y AXIS */ 
    vector<double> x;
    vector<double> y;
    for(i=int(-dim/2); i<int(dim/2); i++){
        for(j=int(-dim/2); j<int(dim/2); j++){
            x.push_back(i*step);
            y.push_back(j*step);
        }
    }

    cout << "Printing Eigenvectors to file" << endl;
    vector<double> efield;
    vector<int> vector_indices;
    for(i=0; i<n; i++){
        if(wr[i] <= k*k*3.5*3.5 && wr[i] > wr[0]){
           vector_indices.push_back(i);
        } 
    }

    /*
    while(vector_indices.empty() == false){
        j = vector_indices.back();
        for(i=j*n; i<(j+1)*n; i++){
            efield.push_back(vr[i]);
        }
        sprintf(vectorfile, "./finemesh/evectors%d.txt", j);
        print_to_file_3d(vectorfile, x, y, efield, n);
        vector_indices.pop_back();
        efield.clear();
    }*/

    /* SAVES FIRST 30 POSSIBLE EVECTORS */
    for(j=0; j<30; j++){
        for(i=j*n; i<(j+1)*n; i++){
            efield.push_back(vr[i]);
        }
        sprintf(vectorfile, "./finemesh/evectors%d.txt", j);
        print_to_file_3d(vectorfile, x, y, efield, n);
        efield.clear();
    }

    /* SAVES EIGENVALUES */
    ofstream myfile("./evalues.txt");
    for(i=0; i<n; i++){
        myfile << wr[i] << "\t" << wi[i] << endl;
    }
    myfile.close();
    
    return 0;
}

void print_to_file_3d(const char* filename, vector<double> &x, vector<double> &y, vector<double> &z, int size){
    /* 
     * Prints data stored in var1 and var2 in two column format in the 
     * file called 'filename', where the length of the arrays var1/2 must
     * be the same.
     * Size is the size of the arrays var1/2 
     * */
    ofstream myfile(filename);
    int i = 0;
    for(i=0; i<size; i++){
        myfile << x[i] << "\t" << y[i] << "\t" << z[i] << endl;
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
                eps[i*dim + j] = n_out*n_out;
            }else if(i>right || j>right){
                eps[i*dim + j] = n_out*n_out;
            }else{
                eps[i*dim + j] = n_in*n_in;
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
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] +  k*k*eps[i%n];
        }
        else if(i ==n*n-1){
            a3[i] = 1/(dy*dy) * 2 * eps[(i-1)%n]/(eps[i%n] + eps[(i-1)%n]);
            a4[i] = 1/(dy*dy);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k*k*eps[i%n];
        }
        else{
            a3[i] = 1/(dy*dy) * 2 * eps[(i-1)%n]/(eps[i%n] + eps[(i-1)%n]);
            a4[i] = 1/(dy*dy) * 2 * eps[(i+1)%n]/(eps[i%n] + eps[(i+1)%n]);
            a5[i] = -2*a1 - 4/(dy*dy) + a3[i] + a4[i] + k*k*eps[i%n];
        }
    }

    
    /* CONTRUCTING MATRIX */
    for(i=0; i < n; i++){
        if( i == 0){
            a[i + i] = a5[i + i];
            a[i*n + i+1] = a2;//4[i*n + i+1];
            a[i*n + i + dim] = a4[i*n + i + dim];//a1
        }
        else if(i == n -1){
            a[i*n + i] = a5[i*n + i];
            a[i*n + i-1] = a1;//3[i*n + i -1];
            a[i*n + i - dim] = a1;//3[i*n + i - dim];
        }else if(i < dim){
            a[i*n + i -1] = a1;//3[i*n + i -1];
            a[i*n + i] = a5[i*n + i];
            a[i*n + i+1] = a2;//4[i*n + i+1];
            a[i*n + i + dim] = a4[i*n + i + dim];//2;            
        }else if(i >= n - dim){
            a[i*n + i - dim] = a3[i*n + i -dim];//1;
            a[i*n + i -1] = a1;//3[i*n + i -1];
            a[i*n + i] = a5[i*n + i];
            a[i*n + i+1] = a2;//4[i*n + i+1];
        }else{
            a[i*n + i - dim] = a3[i*n + i - dim];//1;
            a[i*n + i -1] = a1;//3[i*n + i -1];
            a[i*n + i] = a5[i*n + i];
            a[i*n + i+1] = a2;//4[i*n + i+1];
            a[i*n + i + dim] = a4[i*n + i + dim];//2; 
        }
    }
}  
