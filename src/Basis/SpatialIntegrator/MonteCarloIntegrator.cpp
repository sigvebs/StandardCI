/* 
 * File:   MonteCarloIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 9, 2013, 7:42 PM
 */

#include "MonteCarloIntegrator.h"

//------------------------------------------------------------------------------

MonteCarloIntegrator::MonteCarloIntegrator() {
}

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

double MonteCarloIntegrator::integrate(int p, int q, int r, int s) {
    double I;
    double x, y;
    I = 0;
    for (int i = 0; i < mcSamples; i++) {
        x = (ran3(&idum) - 0.5) * L;
        y = (ran3(&idum) - 0.5) * L;
        I += exp(-w * (x*x + y*y)) * hermitePolynomial(p, sqrtW * x) * hermitePolynomial(q, sqrtW * y) * hermitePolynomial(r, sqrtW * x) * hermitePolynomial(s, sqrtW * y) / (sqrt(pow(x - y, 2) + aa));
    }
    I *= w / (PI * sqrt(pow(2, p + q + r + s) * factorial(p) * factorial(q) * factorial(r) * factorial(s)));
    I /= mcSamples;
    I *= L*L;

#if DEBUG
    cout << "p = " << p << " q = " << q << " r = " << " s = " << s << endl;
    cout << "MonteCarloIntegrator \t= " << I << "\t L = " << L << "\t MC samples = " << mcSamples << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------

MonteCarloIntegrator::MonteCarloIntegrator(const MonteCarloIntegrator& orig) {
}

//------------------------------------------------------------------------------

MonteCarloIntegrator::~MonteCarloIntegrator() {
}
//------------------------------------------------------------------------------
