#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <vector>

#define pi 3.14159265358979323846 

using namespace std;

extern "C" void dstev_(char* job, int* N, double* D, double* OFFD, double* EV, int* VDIM, double* WORK, int* INFO);

void print_to_file(const char* , vector<double> & , vector<double> &, int);
double n_sq(int, int, int, double, double);
double gaussian_sq(double);
double parabola_sq(double);


int main(int argc, char* argv[]){

    int i = 0;
    int mode = 0;
    char modefile[30];
    char indicesfile[30];
    char store_evalues[5];
    
    // parameters for lapack
    char job = 'V';
    int N = 1600, info;

    // systems parameters
    double width = 1.6e-5;
    double scale = 1.0;
    double k = 2*pi/1.55e-6;
    double n1_sq = pow(3.0, 2), n2_sq = pow(3.2, 2);
    double step = width/N;
    double start = -width/2;
    double factor = 1/(k*k*step*step);
    double ratio = 0.5;                         // defines ratio of domain to be given to waveguide (fixed refractive index only)

    if(argc == 2){
        ratio = atof(argv[1]);
        if(ratio > 1 || ratio < 0.1){
            ratio = 0.5;
        }
    }
    cout << "Ratio = " << ratio << endl;

    // create and initialize vector arrays for lapack routine
    vector<double> d, offd, work;
    vector< double > ev(N*N);

    // set up left and right limits based on ratio above
    int left, right;
    left = N*(1-ratio)/2;
    right = N*(1+ ratio)/2;
    
    // set up diagonal elements
    d.push_back(parabola_sq(start) - 2*factor);              // left boundary 
    for(i=1; i<N-1; i++){
        d.push_back(parabola_sq(start + i*step) - 2*factor);
    }
    d.push_back(parabola_sq(-1*start) -2*factor);              // right boundary
   
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

        sprintf(modefile, "./TE/narrow_para_E%d.txt", mode);
        // creates subvector containing one eigenvector
        vector<double> efield(ev.begin() + N*(N-mode-1), ev.begin() + N*(N-mode));
        
        // create vector to represent horizontal plotting axis
        vector<double> distance;
        for(i=int(-N/2); i<int(N/2); i++){
            distance.push_back(i*step);
        }

        // function to print arrays to file for plotting
        print_to_file(modefile, distance, efield, N);
    }

    cout << "Print E-values to file? (y/n)  ";
    cin >> store_evalues;

    if(store_evalues[0] == 'y'){
        vector<double> empty;
        sprintf(indicesfile, "./TE/narrow_parabola_neff.txt");
        print_to_file(indicesfile, d, empty, 20); 
    }

    return 0;
}


void print_to_file(const char* filename, vector<double> &var1, vector<double> &var2, int size){
/* 
Prints data stored in var1 and var2 in two column format in the 
file called 'filename', where the length of the arrays var1/2 must
be the same.
Size is the size of the arrays var1/2 
*/
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

double gaussian_sq(double x){
/*
 Returns squared gaussian which is centred at N/2
 */
    return 9.0*exp(-pow(x, 2)*5e10)*exp(-pow(x, 2)*5e10);
    //return 9.0*exp(-pow(x-0.5*N, 2)*5e10)*exp(-pow(x-0.5*N, 2)*5e10);
}

double parabola_sq(double x){
/*
 Returns squared inverse parabola centred at N/2
 Rather than go negative, this function is discontinuous at +/- 2e-6
 where it goes to 1.
 */
    if(x < -2e-6 || x>=2e-6){
        return pow(3.0 - 5e11*x*x, 2);
    }
    else{
        return 1;
    }
}
