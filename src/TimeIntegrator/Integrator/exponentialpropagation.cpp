#include "exponentialpropagation.h"

//------------------------------------------------------------------------------
ExponentialPropagation::ExponentialPropagation(Config *cfg, HamiltonMatrix *H, cx_vec C): TimeIntegrator(cfg, H, C)
{
    n = H->getDim();
}
//------------------------------------------------------------------------------
void ExponentialPropagation::stepForward()
{
    cx_vec CPrev = C;

    H1 = cx_mat(zeros(n,n), -H->evaluate(t+dt)*dt);
    C = exp(H1)*exp(-H1)*CPrev;

    t +=dt;
}
//------------------------------------------------------------------------------
