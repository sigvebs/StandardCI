/* 
 * File:   Integrator.h
 * Author: sigve
 *
 * Created on December 30, 2012, 3:49 PM
 */

#ifndef INTEGRATOR_H
#define	INTEGRATOR_H

#include "GaussLagandre/gauss_legendre.h"
#include "../Orbital/Orbital.h"
#include <math.h>
#include <cstdlib>

using namespace std;

class Integrator {
public:
    Integrator();
    Integrator(const Integrator& orig);
    virtual ~Integrator();
    bool integrate();
    double k(double x, void* data);
private:

};

#endif	/* INTEGRATOR_H */

