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
void strip_loaded_waveguide(vector<double> &, int, double, double, double, double, double, double, double);
void fd_matrix(vector<double> &, vector<double> &, int, double, double, double, char, char);
void multipliers(char , char , vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, double , double , double , vector<double> &);

int main(){
    
    /* PROGRAM PARAMETERS */
    char mode = 'E', field = 'E';
    int dim = 40, i = 0, j=0;
    int n = dim*dim;
    double k = pi/1.55e-6;
    double width = 1.6e-4;
    double step = width/dim;;
    vector<double> eps(n, 0);
    vector<double> test(n, 0);
    char evectorfile[50];
    char evaluefile[50];
    char guide_dir[20];
    char guidetype;

    /* LAPACK PARAMETERS */
    char jobvl = 'N';
    char jobvr = 'V';
    int N = n;
    int lda = n;
    int ldvl = 1;
    int ldvr = n;
    int lwork = 34*n;          // chosen optimally after test runs
    int info = 1;
    vector<double> wr(n);
    vector<double> wi(n);
    vector<double> a(lda*n, 0);
    vector<double> vl(ldvl*n, 0);
    vector<double> vr(ldvr*n, 0);
    vector<double> work(lwork, 0);

    /* WAVEGUIDE PROPERTIES */
    double ratio = 1/9.0;
    double n_out = 3.16, n_in = 3.5;

    cout << "Construct strip waveguide or core waveguide (s/c)? ";
    cin >> guidetype;
    if(guidetype != 's' && guidetype != 'c'){
        cout << "Invalid waveguide option. Defaulted to strip waveguide" << endl;
        guidetype = 's';
    }

    cout << "Constructing FD Matrix" << endl;    
    if(guidetype == 'c'){
        square_guide_setup(eps, dim, ratio, n_out, n_in);
    }else if(guidetype == 's'){
        strip_loaded_waveguide(eps, dim, width, 2.5e-5, width, 1.5e-5 , 0.5e-5, n_out, n_in);
        cout << "Strip loaded :D" << endl;
    }else{
        cout << "Something went wrong with the waveguide" << endl;
        return 1;
    }
    fd_matrix(a, eps, dim, step, step, k, mode, field);

/*    
    for(i=0; i<8*8; i++){
        if(i%8 == 0){cout << endl;}       
        cout << test[i] << "\t";
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

    if(info == 0){
        cout << "\nSolution Found" << endl;
        cout << "Time Taken: " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
        cout << "Optimal LWORK  = " << work[0] << endl;
    }else{
        cout << "LAPACK Failed\nFORTRAN routine exited with INFO = "<< info << endl;
        cout << "Time Taken: " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
        return 1;
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


    if(guidetype == 's'){
        sprintf(guide_dir, "strip");
    }else{
        sprintf(guide_dir, "core");
    }

    /* SAVES FIRST 30 POSSIBLE EVECTORS */
    for(j=0; j<30; j++){
        for(i=j*n; i<(j+1)*n; i++){
            efield.push_back(vr[i]);
        }
        sprintf(evectorfile, "./coarse_mesh/%s/T%c/%c%d.txt", guide_dir, mode, field, j);
        print_to_file_3d(evectorfile, x, y, efield, n);
        efield.clear();
    }

    /* SAVES EIGENVALUES */
    sprintf(evaluefile, "./coarse_mesh/%s/T%c/evalues.txt", guide_dir, mode);
    ofstream myfile(evaluefile);
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

void strip_loaded_waveguide(vector<double> &eps, int dim, double total_w, double strip_w, double total_h, double strip_h, double strip_base_h, double n_out, double n_in){
    
    int left = dim*(total_w - strip_w)/(2*total_w);
    int right = dim*(total_w + strip_w)/(2*total_w);
    int top = dim*(total_h - strip_h - strip_base_h)/(2*total_h);
    int mid = dim*(total_h + strip_h - strip_base_h)/(2*total_h);
    int bottom = dim*(total_h + strip_h + strip_base_h)/(2*total_h);
    int i,j ;
    for(i=0; i< dim; i++){
        for(j = 0; j< dim; j++){
            if(i< top){
                eps[i*dim + j] = 1;
            }else if(i >= top && i < mid){
                if(j < left || j >= right){eps[i*dim + j] = 1;}
                else{eps[i*dim + j] = n_in;}
            }else if(i >= mid && i < bottom){
                eps[i*dim + j] = n_in;
            }else{
                eps[i*dim + j] = n_out;
            }
        }
    }        
}


void fd_matrix(vector<double> &a, vector<double> &eps, int dim, double dx, double dy, double k, char mode, char field){
/* 
 * CONSTRUCTS FINITE DIFFERENCE MATRIX
 *
 * IN:
 *      a -> vector of dimension dim^4 filled with just zeros (represents FD matrix)
 *      eps -> matrix containing epsilon_r at each coordinate (x, y)
 *      dim -> dimension in of of 1D problem
 *      dx -> step in x direction
 *      dy -> step in y direction
 *      k -> wave vector
 *      mode (E/M) -> represents whether TE or TM mode is wanted
 *      field (E/H) -> represents whether you are interested in the E field or the H field
 *
 * OUT:
 *      a -> vector now filled with elements of FD matrix to be passed into the lapack routine
 *           stored as a vector in order to be compatible with the LAPACK routine
 * ______________________________________________________________________________________________________
 * NOTE:
 * a is a vector of length n*n that represents a n*n matrix
 * to convert between the two storage methods we use the the node number r = i*n + j for point (i, j)
 *
 * This ratio means that the returned vector will have the form --> E = (E00, E01, E02, ..., E0n, E10, E12, ..., ..., Enn)
 * where the first and second indices represent x and y respectively
 */
    int n = dim * dim;
    int i = 0, j=0;
    vector<double> a_up(n, 0);
    vector<double> a_down(n, 0);
    vector<double> a_left(n, 0);
    vector<double> a_right(n, 0);
    vector<double> a_mid(n, 0);

    /* INITIALIZING MULTIPLIERS */
    multipliers(mode, field, a_left, a_right, a_up, a_down, a_mid, dx, dy, k, eps);

    
    /*                          CONTRUCTING MATRIX 
     * 
     * Here we use the modulus operator because the system matrix has dim*dim = n elements
     * but the Finite Difference matrix has n*n elements
     */
    for(i=0; i < n; i++){
        if( i == 0){
            a[i + i] = a_mid[(i + i)%n];
            a[i*n + i+1] = a_up[(i*n + i+1)%n];
            a[i*n + i + dim] = a_right[(i*n + i + dim)%n];
        }
        else if(i == n -1){
            a[i*n + i] = a_mid[(i*n + i)%n];
            a[i*n + i-1] = a_down[(i*n + i-1)%n];
            a[i*n + i - dim] = a_left[(i*n + i -dim)%n];
        }else if(i < dim){
            a[i*n + i -1] = a_down[(i*n + i-1)%n];
            a[i*n + i] = a_mid[(i*n + i)%n];
            a[i*n + i+1] = a_up[(i*n + i+1)%n];
            a[i*n + i + dim] = a_right[(i*n + i + dim)%n];           
        }else if(i >= n - dim){
            a[i*n + i - dim] = a_left[(i*n + i -dim)%n];
            a[i*n + i -1] = a_down[(i*n + i -1)%n];
            a[i*n + i] = a_mid[(i*n + i)%n];
            a[i*n + i+1] = a_up[(i*n + i+1)%n];        
        }else{
            a[i*n + i - dim] = a_left[(i*n + i - dim)%n];
            a[i*n + i -1] = a_down[(i*n + i -1)%n];
            a[i*n + i] = a_mid[(i*n + i)%n];
            a[i*n + i+1] = a_up[(i*n + i+1)%n];
            a[i*n + i + dim] = a_right[(i*n + i + dim)%n]; 
        }
    }
}



void multipliers(char mode_, char field_, vector<double> &a_left, vector<double> &a_right, vector<double> &a_up, vector<double> &a_down, vector<double> &a_mid, double dx, double dy, double k, vector<double> &eps){
/*
 * CALCULATES MULTIPLIERS FOR USE IN FD MATRIX
 *
 * IN:
 *      mode_ (E/M) -> mode you are interested in TM/TE
 *      field_ (E/H) -> field you are interested in 
 *      a_X -> vectors used to store multipliers
 *      dx -> step in x direction
 *      dy -> step in y direction
 *      k -> wave vector
 *      eps -> vector holding relative perimitvities at each point
 *
 * MODIFIED:
 *      a_X -> now contains multipliers to be used when constructing the FD mtrix
 */ 
    
    /* INITIAL CHECK TO ENSURE MULTIPLIERS HAVE A DEFAULT VALUE */
    char mode = mode_;
    char field = field_;
    if(mode != 'E' && mode != 'M'){cout << "---INVALID MODE---" <<endl; mode = 'E';}
    if(field != 'E' && field != 'H'){cout << "---INVALID FIELD---" <<endl; field = 'E';}
    cout << "Mode = " << mode << endl;
    cout << "Field = " << field << endl; 
    
    int n, dim;
    int i;
    n = eps.size();
    dim = sqrt(n);
    if(mode == 'E' && field == 'E'){
        for(i=0; i<n; i++){
            a_up[i] = 1/(dy*dy);
            a_down[i] = 1/(dy*dy);
            if(i%n == 0){
                a_left[i] = 1/(dx*dx);
                a_right[i] = 1/(dx*dx) * 2 * eps[i+1]/(eps[i] + eps[i+1]);
            }
            else if(i%n == 1){
                a_left[i] = 1/(dx*dx) * 2 * eps[i-1]/(eps[i] + eps[i-1]);
                a_right[i] = 1/(dx*dx);
            }
            else{
                a_left[i] = 1/(dx*dx) * 2 * eps[i-1]/(eps[i] + eps[i-1]);
                a_right[i] = 1/(dx*dx) * 2 * eps[i+1]/(eps[i] + eps[i+1]);
            }
            a_mid[i] = -2*a_up[i] - 4/(dx*dx) + a_left[i] + a_right[i] + k*k*eps[i];
        }
    }
    else if(mode == 'M' && field == 'E'){    
        for(i=0; i<n; i++){
            a_right[i] = 1/(dx*dx);
            a_left[i] = 1/(dx*dx);
            if(i < dim){
                a_up[i] = 1/(dy*dy);
                a_down[i] = 1/(dy*dy) * 2 * eps[i+dim]/(eps[i] + eps[i+dim]);
            }
            else if(i >= n-dim){
                a_up[i]  = 1/(dy*dy) * 2 * eps[i-dim]/(eps[i] + eps[i-dim]);
                a_down[i] = 1/(dy*dy);
            }
            else{
                a_up[i]  = 1/(dy*dy) * 2 * eps[i-dim]/(eps[i] + eps[i-dim]);
                a_down[i] = 1/(dy*dy) * 2 * eps[i+dim]/(eps[i] + eps[i+dim]);
            }
            a_mid[i] = -a_left[i] - a_right[i] - 4/(dy*dy) + a_up[i] + a_down[i] +  k*k*eps[i];
        }
    }else if(mode == 'E' && field == 'H'){
        for(i=0; i<n; i++){
            a_up[i] = 1/(dy*dy);
            a_down[i] = 1/(dy*dy);
            if(i%n == 0){
                a_left[i] = 1/(dx*dx);
                a_right[i] = 1/(dx*dx) * 2 * eps[i]/(eps[i] + eps[i+1]);
            }
            else if(i%n == 1){
                a_left[i] = 1/(dx*dx) * 2 * eps[i]/(eps[i] + eps[i-1]);
                a_right[i] = 1/(dx*dx);
            }
            else{
                a_left[i] = 1/(dx*dx) * 2 * eps[i]/(eps[i] + eps[i-1]);
                a_right[i] = 1/(dx*dx) * 2 * eps[i]/(eps[i] + eps[i+1]);
            }
            a_mid[i] =  - a_left[i] - a_right[i] - a_up[i] - a_down[i] + k*k*eps[i];
        }
    }else{
        for(i=0; i<n; i++){
            a_right[i] = 1/(dx*dx);
            a_left[i] = 1/(dx*dx);
            if(i < dim){
                a_up[i] = 1/(dy*dy);
                a_down[i] = 1/(dy*dy) * 2 * eps[i]/(eps[i] + eps[i+dim]);
            }
            else if(i >= n-dim){
                a_up[i]  = 1/(dy*dy) * 2 * eps[i]/(eps[i] + eps[i-dim]);
                a_down[i] = 1/(dy*dy);
            }
            else{
                a_up[i]  = 1/(dy*dy) * 2 * eps[i]/(eps[i] + eps[i-dim]);
                a_down[i] = 1/(dy*dy) * 2 * eps[i]/(eps[i] + eps[i+dim]);
            }
            a_mid[i] = -a_left[i] - a_right[i] - a_up[i] - a_down[i] +  k*k*eps[i];
        }
    }
}
