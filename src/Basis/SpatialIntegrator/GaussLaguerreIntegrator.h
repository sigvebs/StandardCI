/* 
 * File:   GaussLaguerreIntegrator.h
 * Author: sigve
 *
 * Created on January 9, 2013, 8:54 PM
 */

#ifndef GAUSSLAGUERREINTEGRATOR_H
#define	GAUSSLAGUERREINTEGRATOR_H

#include "SpatialIntegrator.h"

class GaussLaguerreIntegrator: public SpatialIntegrator {
public:
    GaussLaguerreIntegrator();
    GaussLaguerreIntegrator(Config *cfg);
    GaussLaguerreIntegrator(const GaussLaguerreIntegrator& orig);
    virtual ~GaussLaguerreIntegrator();
    
    virtual double integrate(int p, int q, int r, int s);
private:
    int N;
    double L;
    double *x;
    double *w1;

};

#endif	/* GAUSSLAGUERREINTEGRATOR_H */

