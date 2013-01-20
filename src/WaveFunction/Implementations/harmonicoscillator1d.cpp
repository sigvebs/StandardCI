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
            * hermitePolynomial(p[1], sqrtW * x[0])
            * hermitePolynomial(q[1], sqrtW * y[0])
            * hermitePolynomial(r[1], sqrtW * x[0])
            * hermitePolynomial(s[1], sqrtW * y[0])
            / (sqrt( pow(x[0] - y[0], 2) + aa));

    return I;
}
//------------------------------------------------------------------------------
double HarmonicOscillator1d::getCoefficient()
{
    return w / (PI
            * sqrt(pow(2, p[1] + q[1] + r[1] + s[1])
            * factorial(p[1])
            * factorial(q[1])
            * factorial(r[1])
            * factorial(s[1])));
}
//------------------------------------------------------------------------------
double HarmonicOscillator1d::getEnergy(vec state)
{
    return w*(state[1] + 0.5);
}
//------------------------------------------------------------------------------
