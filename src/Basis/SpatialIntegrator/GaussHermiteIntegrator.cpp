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
        A = w1[i]
                * hermitePolynomial(p[2], x[i])
                * hermitePolynomial(r[2], x[i]);
        for (int j = 0; j < N; j++) {
            I += A*w1[j]
                    * hermitePolynomial(q[2], x[j])
                    * hermitePolynomial(s[2], x[j])
                    / (sqrt(pow(x[i] - x[j], 2) / w + aa));
        }
    }
    I *=  1/(PI*sqrt(pow(2, p[2] + q[2] + r[2] + s[2]) * factorial(p[2]) * factorial(q[2]) * factorial(r[2]) * factorial(s[2])));

#if DEBUG
    cout << "GH integration = " << I << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------
