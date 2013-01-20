#ifndef HARMONICOSCILLATOR_H
#define HARMONICOSCILLATOR_H

#include "../wavefunction.h"

class HarmonicOscillator: public WaveFunction
{
public:
    HarmonicOscillator(Config *cfg);
    virtual double evaluate(const vec &x, const vec &y) = 0;
    virtual double getCoefficient() = 0;
    virtual double getEnergy(vec state) = 0;
protected:
    double hermitePolynomial(const int degree, const double x);
    int dim;
    double w;
    double sqrtW;
    double aa;
};

#endif // HARMONICOSCILLATOR_H
