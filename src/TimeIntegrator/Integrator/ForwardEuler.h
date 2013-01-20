/* 
 * File:   ForwardEuler.h
 * Author: sigve
 *
 * Created on January 13, 2013, 10:15 PM
 */

#ifndef FORWARDEULER_H
#define	FORWARDEULER_H

#include "../TimeIntegrator.h"

class ForwardEuler: public TimeIntegrator {
public:
    ForwardEuler(Config *cfg, HamiltonMatrix *H, cx_vec C);
    virtual void stepForward();
};

#endif	/* FORWARDEULER_H */

