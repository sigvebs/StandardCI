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

using namespace libconfig;
using namespace std;
using namespace arma;
#include "../../includes/lib.h"


class SpatialIntegrator {
public:
    SpatialIntegrator();
    SpatialIntegrator(Config *cfg);
    SpatialIntegrator(const SpatialIntegrator& orig);
    virtual ~SpatialIntegrator();
    
    virtual double integrate(int p, int q, int r, int s) = 0;

protected:
    double w, sqrtW, a, aa;
    
    double hermitePolynomial(const int degree, const double x);
};

#endif	/* SPATIALINTEGRATOR_H */

