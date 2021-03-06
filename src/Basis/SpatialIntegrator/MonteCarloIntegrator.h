/* 
 * File:   MonteCarloIntegrator.h
 * Author: sigve
 *
 * Created on January 9, 2013, 7:42 PM
 */

#ifndef MONTECARLOINTEGRATOR_H
#define	MONTECARLOINTEGRATOR_H

#include "SpatialIntegrator.h"
#include <math.h>

class MonteCarloIntegrator : public SpatialIntegrator {
public:
    MonteCarloIntegrator(Config *cfg);
    MonteCarloIntegrator(const MonteCarloIntegrator& orig);

    virtual double integrate(const vec &p, const vec &q, const vec &r, const vec &s);
private:
    int mcSamples;
    long idum;
    double L;
};

#endif	/* MONTECARLOINTEGRATOR_H */

