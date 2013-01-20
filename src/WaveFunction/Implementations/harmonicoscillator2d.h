#ifndef HARMONICOSCILLATOR2D_H
#define HARMONICOSCILLATOR2D_H

#include "harmonicoscillator.h"

class HarmonicOscillator2d:public HarmonicOscillator
{
public:
    HarmonicOscillator2d(Config *cfg);
    virtual double evaluate(const vec &r1, const vec &r2);
    virtual double getCoefficient();
    virtual double getEnergy(vec state);
};

#endif // HARMONICOSCILLATOR2D_H
