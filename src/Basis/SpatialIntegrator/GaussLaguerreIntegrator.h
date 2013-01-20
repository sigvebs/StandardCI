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
    GaussLaguerreIntegrator(Config *cfg, WaveFunction *wf);
    GaussLaguerreIntegrator(const GaussLaguerreIntegrator& orig);
    virtual ~GaussLaguerreIntegrator();

    virtual double integrate(const vec &p, const vec &q, const vec &r, const vec &s);
private:
    vec x_, y_;
    int N;
    double L;
    double *x;
    double *w1;
};

#endif	/* GAUSSLAGUERREINTEGRATOR_H */

