#include "harmonicoscillator1d.h"

//------------------------------------------------------------------------------
HarmonicOscillator1d::HarmonicOscillator1d(Config *cfg):HarmonicOscillator(cfg)
{
}
//------------------------------------------------------------------------------
double HarmonicOscillator1d::evaluate(const vec &x, const vec &y)
{
    double I;
    I = exp(-w * (x[0] * x[0] + y[0] * y[0]))
            * hermitePolynomial(p[2], sqrtW * x[0])
            * hermitePolynomial(q[2], sqrtW * y[0])
            * hermitePolynomial(r[2], sqrtW * x[0])
            * hermitePolynomial(s[2], sqrtW * y[0])
            / (sqrt( pow(x[0] - y[0], 2) + aa));

    return I;
}
//------------------------------------------------------------------------------
double HarmonicOscillator1d::getCoefficient()
{
    return w / (PI
            * sqrt(pow(2, p[2] + q[2] + r[2] + s[2])
            * factorial(p[2])
            * factorial(q[2])
            * factorial(r[2])
            * factorial(s[2])));
}
//------------------------------------------------------------------------------
double HarmonicOscillator1d::getEnergy(vec state)
{
    return w*(state[2] + 0.5);
}
//------------------------------------------------------------------------------
