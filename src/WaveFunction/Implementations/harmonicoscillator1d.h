#ifndef HARMONICOSCILLATOR1D_H
#define HARMONICOSCILLATOR1D_H

#include "harmonicoscillator.h"

class HarmonicOscillator1d: public HarmonicOscillator
{
public:
    HarmonicOscillator1d(Config *cfg);
    virtual double evaluate(const vec &x, const vec &y);
    virtual double getCoefficient();
    virtual double getEnergy(vec state);
};

#endif // HARMONICOSCILLATOR1D_H
