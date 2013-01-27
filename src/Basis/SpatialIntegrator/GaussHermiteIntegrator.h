/* 
 * File:   GaussHermiteIntegrator.h
 * Author: sigve
 *
 * Created on January 11, 2013, 4:21 PM
 */

#ifndef GAUSSHERMITEINTEGRATOR_H
#define	GAUSSHERMITEINTEGRATOR_H

#include "SpatialIntegrator.h"

class GaussHermiteIntegrator : public SpatialIntegrator {
public:
    GaussHermiteIntegrator(Config* cfg);
    GaussHermiteIntegrator(const GaussHermiteIntegrator& orig);
    virtual double integrate(const vec &p, const vec &q, const vec &r, const vec &s);
private:
    int N;
    double *x;
    double *w1;
};

#endif	/* GAUSSHERMITEINTEGRATOR_H */

