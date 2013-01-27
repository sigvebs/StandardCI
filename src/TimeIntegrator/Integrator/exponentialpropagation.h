#ifndef EXPONENTIALPROPAGATION_H
#define EXPONENTIALPROPAGATION_H

#include "../TimeIntegrator.h"

class ExponentialPropagation: public TimeIntegrator
{
public:
    ExponentialPropagation(Config *cfg, HamiltonMatrix *H, cx_vec C);
    virtual void stepForward();
private:
    cx_mat H1;
    int n;
};

#endif // EXPONENTIALPROPAGATION_H
