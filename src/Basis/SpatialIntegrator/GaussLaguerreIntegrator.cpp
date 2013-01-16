/* 
 * File:   GaussLaguerreIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 9, 2013, 8:54 PM
 */

#include "GaussLaguerreIntegrator.h"

//------------------------------------------------------------------------------

GaussLaguerreIntegrator::GaussLaguerreIntegrator() {
}
//------------------------------------------------------------------------------

GaussLaguerreIntegrator::GaussLaguerreIntegrator(Config *cfg) : SpatialIntegrator(cfg) {

    try {
        cfg->lookupValue("spatialIntegration.GaussLaguerre.samples", N);
        cfg->lookupValue("spatialIntegration.L", L);
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Config *cfg)::Error reading from 'systemSettings' object setting." << endl;
    }

#if DEBUG
    cout << "GaussLaguerreIntegrator(Setting* systemSettings)" << endl;
    cout << "N = " << N << endl;
    cout << "L = " << L << endl;
#endif

    x = new double[N];
    w1 = new double[N];
    gauleg(-L, L, x, w1, N);
}

//------------------------------------------------------------------------------

double GaussLaguerreIntegrator::integrate(int p, int q, int r, int s) {
    double I = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            I += w1[i] * w1[j]
                    * exp(-w * (x[i] * x[i] + x[j] * x[j]))
                    * hermitePolynomial(p, sqrtW * x[i]) * hermitePolynomial(q, sqrtW * x[j])
                    * hermitePolynomial(r, sqrtW * x[i]) * hermitePolynomial(s, sqrtW * x[j])
                    / (sqrt(pow(x[i] - x[j], 2) + aa));
        }
    }
    I *= w / (PI * sqrt(pow(2, p + q + r + s) * factorial(p) * factorial(q) * factorial(r) * factorial(s)));

#if DEBUG
    cout << "GL integration = " << I << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------

GaussLaguerreIntegrator::GaussLaguerreIntegrator(const GaussLaguerreIntegrator& orig) {
}

//------------------------------------------------------------------------------

GaussLaguerreIntegrator::~GaussLaguerreIntegrator() {
    delete x;
    delete w1;
}

//------------------------------------------------------------------------------
