/* 
 * File:   CrankNicolson.h
 * Author: sigve
 *
 * Created on January 14, 2013, 10:06 AM
 */

#ifndef CRANKNICOLSON_H
#define	CRANKNICOLSON_H

#include "../TimeIntegrator.h"

class CrankNicolson: public TimeIntegrator {
public:
    CrankNicolson();
    CrankNicolson(Config *cfg, HamiltonMatrix *H, cx_vec C);
    CrankNicolson(const CrankNicolson& orig);
    virtual ~CrankNicolson();
    
    virtual void stepForward();
private:
    cx_mat H1;
    cx_mat H2;
    cx_mat I;
    int n;

};

#endif	/* CRANKNICOLSON_H */

