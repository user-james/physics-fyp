#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <vector>

#define pi 3.14159265358979323846 

using namespace std;

extern "C" void dstev_(char* job, int* N, double* D, double* OFFD, double* EV, int* VDIM, double* WORK, int* INFO);

void print_to_file(const char* , double* , double*, int);
double n_sq(int, int, int, double, double);
double gaussian(int, int);


int main(){

    int i = 0;
    int mode = 0;
    char filename[30];
    
    // parameters for lapack
    char job = 'V';
    int N = 1600, info;

    // systems parameters
    double a = 2e10-6;
    double k = pi/1.55e-6;
    double n1_sq = pow(3.0, 2), n2_sq = pow(3.1, 2);
    double step = 1e-7;
    double factor = 1/(k*k*step*step);
    double ratio = 0.05;                         // defines ratio of domain to be given to waveguide (fixed refractive index only)

    // create and initialize vector arrays for lapack routine
    vector<double> d, offd, work;
    vector< double > ev(N*N);

    // set up left and right limits based on ratio above
    int left, right;
    left = N*(1-ratio)/2;
    right = N*(1+ ratio)/2;
    
    // set up diagonal elements
    d.push_back(gaussian(0, N) - 2*factor);              // left boundary 
    for(i=1; i<N-1; i++){
        d.push_back(gaussian(i, N) - 2*factor);
    }
    d.push_back(gaussian(N, N) - 2*factor);              // right boundary
   
    cout << "Gaussian(N) = " << gaussian(0.5*N, N) << endl; 
    // set up off diagonal elements
    for(i=0; i<N-1; i++){
        offd.push_back(factor);
    }
   
    // initialize other arrays
    for(i=0; i<2*N-2; i++){
        work.push_back(0.0);
    }

    cout << "Solving for Eigenvalues/vectors, may take a while\n"; 

    // run FORTRAN routine to solve for eigenvalues/vectors
    dstev_(&job, &N, &d[0], &offd[0], &ev[0], &N, &*work.begin(), &info);
   
    cout << "Solution found\n";

    while(true){
        cout << "\nSelect mode number you are interested in: ";
        cin >> mode;

        if(mode < 0){break;}

        sprintf(filename, "./data/TE%d.txt", mode);
        // creates subvector containing one eigenvector
        vector<double> efield(ev.begin() + N*(N-mode-1), ev.begin() + N*(N-mode));
        
        // create vector to represent horizontal plotting axis
        vector<double> distance;
        for(i=int(-N/2); i<int(N/2); i++){
            distance.push_back(i*step);
        }

        // function to print arrays to file for plotting
        print_to_file(filename, &distance[0], &efield[0], N);
    }

    return 0;
}


void print_to_file(const char* filename, double *var1, double *var2, int size){
/* 
Prints data stored in var1 and var2 in two column format in the 
file called 'filename', where the length of the arrays var1/2 must
be the same.
Size is the size of the arrays var1/2 
*/
    ofstream myfile(filename);
    int i = 0;

    for(i=0; i<size; i++){
        myfile << var1[i] << "\t" << var2[i] << endl;
    }

    myfile.close();

}

double n_sq(int x, int left, int right, double n1_sq, double n2_sq){
/*
Returns refractive index function
*/
    if(x < left || x>= right){
        return n1_sq;
    }
    else{
        return n2_sq;
    }

}

double gaussian(int x, int N){
/*
 Returns squared gaussian which is centred at N/2
 */
    return 9.0*exp(-pow(x-0.5*N, 2)*5e10)*exp(-pow(x-0.5*N, 2)*5e10);
}
