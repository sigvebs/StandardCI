/* 
 * File:   Integrator.cpp
 * Author: sigve
 * 
 * Created on December 30, 2012, 3:49 PM
 */

#include <iostream>

#include "Integrator.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#ifndef FABS
#define FABS(a) ((a)>=0?(a):-(a))
#endif

//------------------------------------------------------------------------------

double f(double x, void* data) {
    return sin(x);
}

//------------------------------------------------------------------------------

double Integrator::k(double x, void* data) {
    return sin(x);
}

//------------------------------------------------------------------------------

Integrator::Integrator() {
}

//------------------------------------------------------------------------------

Integrator::Integrator(const Integrator& orig) {
}

//------------------------------------------------------------------------------

Integrator::~Integrator() {
}

//------------------------------------------------------------------------------

bool Integrator::integrate() {
    /* numerical approximation of integral */
    double approx;

    /* true value of int(sin(x), x=0..Pi) = 2.0*/
    double exact = 2.0;

    /* approximation error */
    double error;

    int i;

    cout << "Numerical Approximation of int(sin(x), x=0..Pi) by Gauss-Legendre Quadrature:\n" << endl;

    for (i = 2; i <= 128; i++) {
        approx = gauss_legendre(i, f, NULL, 0, PI);
        error = approx - exact;
        //        cout << "n = " << i << " , approx = " << approx << endl;
        printf("n = %4d: error = %.15g\n", i, FABS(error));
    }
}

//------------------------------------------------------------------------------