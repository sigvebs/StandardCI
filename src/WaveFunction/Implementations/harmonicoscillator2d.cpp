#include "harmonicoscillator2d.h"

//------------------------------------------------------------------------------
HarmonicOscillator2d::HarmonicOscillator2d(Config *cfg):HarmonicOscillator(cfg)
{
    aa = 0.25*0.25;
}
//------------------------------------------------------------------------------
double HarmonicOscillator2d::evaluate(const vec &r1, const vec &r2)
{
    double I;

    I = exp(-w * (dot(r1,r1) + dot(r2,r2)));

    I *=      hermitePolynomial(p[2], sqrtW * r1[0])
            * hermitePolynomial(p[3], sqrtW * r1[1])

            * hermitePolynomial(q[2], sqrtW * r2[0])
            * hermitePolynomial(q[3], sqrtW * r2[1])

            * hermitePolynomial(r[2], sqrtW * r1[0])
            * hermitePolynomial(r[3], sqrtW * r1[1])

            * hermitePolynomial(s[2], sqrtW * r2[0])
            * hermitePolynomial(s[3], sqrtW * r2[1])

            / sqrt(dot(r1 - r2,r1 - r2) + aa);

    return I;
}
//------------------------------------------------------------------------------
double HarmonicOscillator2d::getCoefficient()
{
    double psi = w / (PI
            * sqrt(pow(2, p[2] + q[2] + r[2] + s[2] + p[3] + q[3] + r[3] + s[3])
            * factorial(p[2]) * factorial(q[2]) * factorial(r[2]) * factorial(s[2])
            * factorial(p[3]) * factorial(q[3]) * factorial(r[3]) * factorial(s[3])));

    return psi;
}
//------------------------------------------------------------------------------
double HarmonicOscillator2d::getEnergy(vec state)
{
//    return w*(state[2] + state[3] + 1);
    return w*(2*state[2] + abs(state[3]) + 1);
}
//------------------------------------------------------------------------------
