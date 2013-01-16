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

double GaussHermiteIntegrator::integrate(int p, int q, int r, int s) {
    double I = 0;
    double A;
    for (int i = 0; i < N; i++) {
        A = w1[i] * hermitePolynomial(p, x[i]) * hermitePolynomial(r, x[i]);
        for (int j = 0; j < N; j++) {
            I += A * w1[j]
                    * hermitePolynomial(q, x[j])
                    * hermitePolynomial(s, x[j])
                    / (sqrt(pow(x[i] - x[j], 2) / w + aa));
        }
    }
    I *= 1 / (PI * sqrt(pow(2, p + q + r + s) * factorial(p) * factorial(q) * factorial(r) * factorial(s)));

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
