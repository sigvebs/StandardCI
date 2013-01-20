/* 
 * File:   GaussHermiteIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 11, 2013, 4:21 PM
 */

#include "GaussHermiteIntegrator.h"

//------------------------------------------------------------------------------

GaussHermiteIntegrator::GaussHermiteIntegrator() {
}

//------------------------------------------------------------------------------

GaussHermiteIntegrator::GaussHermiteIntegrator(Config *cfg) : SpatialIntegrator(cfg) {

    try {
        cfg->lookupValue("spatialIntegration.GaussHermite.samples", N);
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Setting* systemSettings)::Error reading from 'systemSettings' object setting." << endl;
    }

#if DEBUG
    cout << "GaussLaguerreIntegrator(Config *cfg)" << endl;
    cout << "N = " << N << endl;
#endif

    x = new double[N];
    w1 = new double[N];

    gauher(x, w1, N);
    
//    for(int i=0; i<N; i++){
//        cout << "x[" << i << "] = " << x[i] << endl;
//    }
//    for(int i=0; i<N; i++){
//        cout << "w[" << i << "] = " << w1[i] << endl;
//    }

}
//------------------------------------------------------------------------------

double GaussHermiteIntegrator::integrate(const vec &p, const vec &q, const vec &r, const vec &s) {
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

//    cout << "GH integration = " << I << endl;
#if DEBUG
    cout << "GH integration = " << I << endl;
#endif
    return I;
}

//------------------------------------------------------------------------------
GaussHermiteIntegrator::GaussHermiteIntegrator(const GaussHermiteIntegrator & orig) {
}
//------------------------------------------------------------------------------
GaussHermiteIntegrator::~GaussHermiteIntegrator() {
}
//------------------------------------------------------------------------------
