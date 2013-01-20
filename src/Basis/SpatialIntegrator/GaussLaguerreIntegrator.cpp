/* 
 * File:   GaussLaguerreIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 9, 2013, 8:54 PM
 */

#include "GaussLaguerreIntegrator.h"

//------------------------------------------------------------------------------
GaussLaguerreIntegrator::GaussLaguerreIntegrator(Config *cfg, WaveFunction *wf) : SpatialIntegrator(cfg, wf)
{
    try {
        cfg->lookupValue("spatialIntegration.GaussLaguerre.samples", N);
        cfg->lookupValue("spatialIntegration.L", L);
        dim = cfg->lookup("systemSettings.dim");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Config *cfg)::Error reading from 'systemSettings' object setting." << endl;
    }

    x_ = zeros(dim);
    y_ = zeros(dim);
    x = new double[N];
    w1 = new double[N];
    gauleg(-L, L, x, w1, N);

#if DEBUG
    cout << "GaussLaguerreIntegrator(Setting* systemSettings)" << endl;
    cout << "N = " << N << endl;
    cout << "L = " << L << endl;
#endif
}
//------------------------------------------------------------------------------
double GaussLaguerreIntegrator::integrate(const vec &p, const vec &q, const vec &r, const vec &s)
{
    double I = 0;
    wf->setQuantumNumber(p,q,r,s);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            x_[0] = x[i];
            y_[0] = x[j];
            I += w1[i] * w1[j] * wf->evaluate(x_, y_);
        }
    }
    I *= wf->getCoefficient();
#if DEBUG
    cout << "GL integration = " << I << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------
GaussLaguerreIntegrator::~GaussLaguerreIntegrator()
{
    delete x;
    delete w1;
}
//------------------------------------------------------------------------------
