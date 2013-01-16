/* 
 * File:   Orbital.h
 * Author: sigve
 *
 * Created on December 28, 2012, 3:24 PM
 */

#ifndef ORBITAL_H
#define	ORBITAL_H
#include <armadillo>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;
using namespace arma;

class Orbital {
public:
    Orbital();
    Orbital(vec quantumNumbers);
    Orbital(const Orbital& orig);
    virtual ~Orbital();
    virtual double evaluate(vec x) = 0;
private:
   vec quantumNumbers;
};

#endif	/* ORBITAL_H */

