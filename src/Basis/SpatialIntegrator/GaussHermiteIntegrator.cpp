/* 
 * File:   GaussHermiteIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 11, 2013, 4:21 PM
 */

#include "GaussHermiteIntegrator.h"


//------------------------------------------------------------------------------
GaussHermiteIntegrator::GaussHermiteIntegrator(Config *cfg) : SpatialIntegrator(cfg)
{
    try {
        cfg->lookupValue("spatialIntegration.GaussHermite.samples", N);
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Setting* systemSettings)::Error reading from 'systemSettings' object setting." << endl;
    }

    x = new double[N];
    w1 = new double[N];

    gauher(x, w1, N);
#if DEBUG
    cout << "GaussLaguerreIntegrator(Config *cfg)" << endl;
    cout << "N = " << N << endl;
#endif
}
//------------------------------------------------------------------------------
double GaussHermiteIntegrator::integrate(const vec &p, const vec &q, const vec &r, const vec &s)
{
    double I = 0;
    double A;
    for (int i = 0; i < N; i++) {
        A = w1[i] * hermitePolynomial(p[1], x[i]) * hermitePolynomial(r[1], x[i]);
        for (int j = 0; j < N; j++) {
            I += A * w1[j]
                    * hermitePolynomial(q[1], x[j])
                    * hermitePolynomial(s[1], x[j])
                    / (sqrt(pow(x[i] - x[j], 2) / w + aa));
        }
    }
    I *= 1 / (PI * sqrt(pow(2, p[1] + q[1] + r[1] + s[1]) * factorial(p[1]) * factorial(q[1]) * factorial(r[1]) * factorial(s[1])));

#if DEBUG
    cout << "GH integration = " << I << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------
