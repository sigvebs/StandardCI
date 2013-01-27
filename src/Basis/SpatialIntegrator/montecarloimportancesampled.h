#ifndef MONTECARLOIMPORTANCESAMPLED_H
#define MONTECARLOIMPORTANCESAMPLED_H

#include "SpatialIntegrator.h"

class MonteCarloImportanceSampled : public SpatialIntegrator
{
public:
    MonteCarloImportanceSampled(Config *cfg);

    virtual double integrate(const vec &p, const vec &q, const vec &r, const vec &s);
private:
    int mcSamples;
    long idum;
    double L;
};

#endif // MONTECARLOIMPORTANCESAMPLED_H
