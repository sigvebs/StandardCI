#include "wavefunction.h"

//------------------------------------------------------------------------------
WaveFunction::WaveFunction(Config *cfg): cfg(cfg)
{
}
//------------------------------------------------------------------------------
void WaveFunction::setQuantumNumber(const vec &p_, const vec &q_, const vec &r_, const vec &s_)
{
    p = p_;
    q = q_;
    r = r_;
    s = s_;
}
//------------------------------------------------------------------------------
