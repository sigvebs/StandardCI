#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <libconfig.h++>
#include "../includes/lib.h"
#include "../headers.h"

using namespace std;
using namespace libconfig;
using namespace arma;

class WaveFunction
{
public:
    WaveFunction(Config *cfg);
    virtual double evaluate(const vec &x, const vec &y) = 0;
    virtual double getCoefficient() = 0;
    virtual double getEnergy(vec state) = 0;
    void setQuantumNumber(const vec &p_, const vec &q_, const vec &r_, const vec &s_);
protected:
    Config *cfg;
    vec p,q,r,s;
};

#endif // WAVEFUNCTION_H
