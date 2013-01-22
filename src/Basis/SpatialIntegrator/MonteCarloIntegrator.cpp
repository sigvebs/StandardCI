/* 
 * File:   MonteCarloIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 9, 2013, 7:42 PM
 */

#include "MonteCarloIntegrator.h"

//------------------------------------------------------------------------------
MonteCarloIntegrator::MonteCarloIntegrator(Config *cfg) : SpatialIntegrator(cfg) {

    try {
        cfg->lookupValue("spatialIntegration.MonteCarlo.samples", mcSamples);
        cfg->lookupValue("spatialIntegration.L", L);
        int dummy;
        cfg->lookupValue("spatialIntegration.idum", dummy);
        idum = dummy;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Setting* systemSettings)::Error reading from 'systemSettings' object setting." << endl;
    }

#if DEBUG
    cout << "MonteCarloIntegrator::MonteCarloIntegrator(Config *cfg)" << endl;
    cout << "mcSamples = " << mcSamples << endl;
    cout << "L = " << L << endl;
    cout << "idum = " << idum << endl;
#endif
}
//------------------------------------------------------------------------------
double MonteCarloIntegrator::integrate(const vec &p, const vec &q, const vec &r, const vec &s) {
    double I;
    double x, y;
    I = 0;
    for (int i = 0; i < mcSamples; i++) {
        x = (ran3(&idum) - 0.5) * L;
        y = (ran3(&idum) - 0.5) * L;
        I += exp(-w * (x*x + y*y))
                * hermitePolynomial(p[2], sqrtW * x)
                * hermitePolynomial(q[2], sqrtW * y)
                * hermitePolynomial(r[2], sqrtW * x)
                * hermitePolynomial(s[2], sqrtW * y)
                / (sqrt(pow(x - y, 2) + aa));
    }
    I *= w / (PI * sqrt(pow(2, p[2] + q[2] + r[2] + s[2]) * factorial(p[2]) * factorial(q[2]) * factorial(r[2]) * factorial(s[2])));
    I /= mcSamples;
    I *= L*L;

#if DEBUG
    cout << "p = " << p[2] << " q = " << q[2] << " r = " << r[2] << " s = " << s[2] << endl;
    cout << "MonteCarloIntegrator \t= " << I << "\t L = " << L << "\t MC samples = " << mcSamples << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------
