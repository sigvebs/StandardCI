/* 
 * File:   SpatialIntegrator.h
 * Author: sigve
 *
 * Created on January 9, 2013, 7:37 PM
 */

#ifndef SPATIALINTEGRATOR_H
#define	SPATIALINTEGRATOR_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <libconfig.h++>

#include "../../headers.h"
#include "../../WaveFunction/wavefunction.h"

using namespace libconfig;
using namespace std;
using namespace arma;
#include "../../includes/lib.h"


class SpatialIntegrator {
public:
    SpatialIntegrator(Config *cfg);
    SpatialIntegrator(Config *cfg, WaveFunction *wf);
    SpatialIntegrator(const SpatialIntegrator& orig);
    virtual double integrate(const vec &p, const vec &q, const vec &r, const vec &s) = 0;
protected:
    WaveFunction *wf;
    int dim;
    double w, sqrtW, a, aa;
    
    // TODO: should not be in this class. Remove when WF has been implemented
    // trough the entire program.
    double hermitePolynomial(const int degree, const double x);
};

#endif	/* SPATIALINTEGRATOR_H */

