/* 
 * File:   MonteCarloIntegrator.h
 * Author: sigve
 *
 * Created on January 9, 2013, 7:42 PM
 */

#ifndef MONTECARLOINTEGRATOR_H
#define	MONTECARLOINTEGRATOR_H

#include "SpatialIntegrator.h"

class MonteCarloIntegrator : public SpatialIntegrator {
public:
    MonteCarloIntegrator();
    MonteCarloIntegrator(Config *cfg);
    MonteCarloIntegrator(const MonteCarloIntegrator& orig);
    virtual ~MonteCarloIntegrator();

    virtual double integrate(int p, int q, int r, int s);
private:
    int mcSamples;
    long idum;
    double L;
};

#endif	/* MONTECARLOINTEGRATOR_H */

