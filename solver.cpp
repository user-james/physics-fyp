#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>

#define PI 3.14159265359

using namespace std;

double bisection(double, double, double (*f)(double, double), double, double);
double solve(double, double (*f)(double, double), double, double, double);
double newton(double, double (*f)(double, double), double (*df)(double, double), double, double);
double func(double, double);
double dfunc(double, double);
int sign(double);
double yfunc(double);


int main(){

    // initial parameters
    double lambda = 1.55e-6, n2 = 3.1, n1 = 3.0, a = 2.0e-6;
    double k = 2*PI/lambda;
    double r = k*a*sqrt(n2*n2 - n1*n1);
    double guess = 1.0, x = 0, y = 0, alpha = 0.5, left, right;
    string fname = "tmvalues";

    //cout << "Enter guess: ";
    //cin >> guess;
    //cout << "Enter alpha: ";
    //cin >> alpha;
    cout << "Left limit: ";
    cin >> left;
    cout << "Right limit: ";
    cin >> right;


    //x = solve(guess, func, 1e-7, alpha, r);
    //x = newton(guess, func, dfunc, 1e-10, r);
    x = bisection(left, right, func, 1e-10, r);
    y = yfunc(x);

    cout << "r = " << r << endl;
    cout << "x = " << x << "\ty = " << y << endl;    
    cout << "sgn(x)= " << sign(x) << endl;
    cout << "func(4) = " << func(4, r) << endl;

    fname = fname + to_string(x) + ".txt";
    ofstream param_file;
    param_file.open(fname);
    param_file << r << "\n" << a << '\n' << x << '\n' << y << endl;
    param_file.close();

    return 0;
}


double solve(double guess, double (*f)(double, double), double err,  double alpha, double r){
// solves for root of f(x) near guess
// with alpha as a parameter of the solver and r is a parameter for the function

    double x = 0, prev = guess, temp = 0;
    int count = 1;
    x = (*f)(guess, r);

   
    while(fabs(prev -x) > err ){
        temp = prev;
        prev = x;
        x = (1 - alpha)*(*f)(temp, r) + alpha*(*f)(prev, r);
        count++;

        // including breaking statement 
        if(count > 1000){
            cout << "Max. iterations reached\n";
            break;
        }
    }
    
    cout << "# iterations: " << count << endl;
    return x;

}

double newton(double guess, double (*f)(double, double), double (*df)(double, double), double err, double r){
// solves using newton-raphson

    double x = 0, prev = guess;
    int count = 1;
    x = prev - (*f)(prev, r)/(*df)(prev, r);

    while(fabs(prev -x) > err ){
        prev = x;
        x = prev - (*f)(prev, r)/(*df)(prev, r);
        count++;
        
        // including breaking statement 
        if(count > 1000){
            cout << "Max. iterations reached\n";
            break;
        }
    }
    
    cout << "# iterations: " << count << endl;
    return x;


}

double bisection(double a, double b, double (*f)(double, double), double err, double r){
//solves using the bisection method
    double c;
    int count = 0;

    while(count < 2000){
        c = (a + b)/2.0;
        if(fabs(a-b) <= err || (*f)(c, r) == 0){
            cout << "# iterations: " << count << endl;
            return c;
        }
        else if(sign((*f)(a, r)) == sign((*f)(c, r))){ a = c;}
        else{b = c;}
        count++;
    }

    cout << "Max iterations reached" << endl;
    return 0;
}


double func(double x, double r){
// function in x whose root needs to be found
// r is a parameter
    return x - r *pow(1 + (1.06778/tan(x))*(1.06778/tan(x)), -0.5);
}

int sign(double x){
// returns the sign of x
    return (x > 0) - (x < 0);
}

double dfunc(double x, double r){

    return 1 - r*cos(x);
}

double yfunc(double x){

    return -x*1.06778/tan(x);
}
