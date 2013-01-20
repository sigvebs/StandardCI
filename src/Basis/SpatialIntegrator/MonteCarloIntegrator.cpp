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
                * hermitePolynomial(p[1], sqrtW * x)
                * hermitePolynomial(q[1], sqrtW * y)
                * hermitePolynomial(r[1], sqrtW * x)
                * hermitePolynomial(s[1], sqrtW * y)
                / (sqrt(pow(x - y, 2) + aa));
    }
    I *= w / (PI * sqrt(pow(2, p[1] + q[1] + r[1] + s[1]) * factorial(p[1]) * factorial(q[1]) * factorial(r[1]) * factorial(s[1])));
    I /= mcSamples;
    I *= L*L;

#if DEBUG
    cout << "p = " << p[1] << " q = " << q[1] << " r = " << r[1] << " s = " << s[1] << endl;
    cout << "MonteCarloIntegrator \t= " << I << "\t L = " << L << "\t MC samples = " << mcSamples << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------
