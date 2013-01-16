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
    ForwardEuler();
    ForwardEuler(Config *cfg, HamiltonMatrix *H, cx_vec C);
    ForwardEuler(const ForwardEuler& orig);
    virtual ~ForwardEuler();
    
    virtual void stepForward();
private:
    

};

#endif	/* FORWARDEULER_H */

