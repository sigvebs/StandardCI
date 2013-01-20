/* 
 * File:   TimeIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 13, 2013, 9:51 PM
 */

#include "TimeIntegrator.h"
//------------------------------------------------------------------------------
TimeIntegrator::TimeIntegrator(Config *cfg, HamiltonMatrix *H, cx_vec C) : H(H), C(C)
{
    try {
        N = cfg->lookup("TimeIntegration.N");
        dt = cfg->lookup("TimeIntegration.dt");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "TimeIntegrator(Config *cfg)::Error reading from 'TimeIntegration' object setting." << endl;
    }
    t = 0;
#if DEBUG
    cout << "TimeIntegrator(Config *cfg)" << endl;
    cout << "N = " << N << endl;
    cout << "dt = " << dt << endl;
#endif
}
//------------------------------------------------------------------------------
cx_vec TimeIntegrator::getCoefficients()
{
    return C;
}
//------------------------------------------------------------------------------
