/* 
 * File:   TimeIntegrator.h
 * Author: sigve
 *
 * Created on January 13, 2013, 9:51 PM
 */

#ifndef TIMEINTEGRATOR_H
#define	TIMEINTEGRATOR_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <libconfig.h++>

#include "../headers.h"
#include "../Hamiltonian/HamiltonMatrix.h"

using namespace libconfig;
using namespace std;
using namespace arma;

class TimeIntegrator {
public:
    TimeIntegrator();
    TimeIntegrator(Config *cfg, HamiltonMatrix *H, cx_vec C);
    TimeIntegrator(const TimeIntegrator& orig);
    virtual ~TimeIntegrator();
    
    virtual void stepForward() = 0;
    
    cx_vec getCoefficients();
protected:
    int N;
    double dt, t;
    HamiltonMatrix *H;
    cx_vec C;
};

#endif	/* TIMEINTEGRATOR_H */

