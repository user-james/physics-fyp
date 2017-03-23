#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <vector>

#define pi 3.14159265358979323846 

using namespace std;

extern "C" void dstev_(char* job, int* N, double* D, double* OFFD, double* EV, int* VDIM, double* WORK, int* INFO);

void square_guide_setup(vector< vector<double> > &, int, double, double, double);


int main(){

    int dim = 4, i = 0, j=0;
    
    
    cout << "Enter size of square matrix: ";
    cin >> dim;
    vector< vector<double> > eps;
    eps.resize(dim, vector<double>(dim, 0));

    double ratio = 1/4444.0;
    double n_out = 1.4, n_in = 1.5;

    square_guide_setup(eps, dim, ratio, n_out, n_in);

    for(j=0; j<dim; j++){
        for(i=0; i<dim; i++){ cout << eps[i][j] << "\t";}
        cout << endl;
    }


    return 0;
}



void square_guide_setup(vector< vector<double> > &eps , int dim, double centre_ratio, double n_out, double n_in){
    
    int left, right, i = 0, j = 0;
    left = dim * (1. - centre_ratio) / 2;
    right = dim  - 1 - left;

    for(i=0; i<dim; i++){
        for(j = 0; j<dim; j++){
            if( i<left || j<left){
                eps[i][j] = n_out;
            }else if(i>right || j>right){
                eps[i][j] = n_out;
            }else{
                eps[i][j] = n_in;
            }            
        }
    }

 
}
