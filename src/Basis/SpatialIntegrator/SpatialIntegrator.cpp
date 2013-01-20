/* 
 * File:   SpatialIntegrator.cpp
 * Author: sigve
 * 
 * Created on January 9, 2013, 7:37 PM
 */

#include "SpatialIntegrator.h"
//------------------------------------------------------------------------------
SpatialIntegrator::SpatialIntegrator(Config *cfg) {
    try {
        cfg->lookupValue("potential.a", a);
        cfg->lookupValue("systemSettings.w", w);
        dim = cfg->lookup("systemSettings.dim");
        aa = a*a;
        sqrtW = sqrt(w);
     } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Setting* systemSettings)::Error reading from 'systemSettings' object setting." << endl;
    }
}
//------------------------------------------------------------------------------
SpatialIntegrator::SpatialIntegrator(Config *cfg, WaveFunction *wf): wf(wf) {
    try {
        cfg->lookupValue("potential.a", a);
        cfg->lookupValue("systemSettings.w", w);
        aa = a*a;
        sqrtW = sqrt(w);
     } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Setting* systemSettings)::Error reading from 'systemSettings' object setting." << endl;
    }
}
//------------------------------------------------------------------------------
double SpatialIntegrator::hermitePolynomial(const int degree, const double x) {
    double hermite_polynomial;

    switch (degree) {
        case 0:
            hermite_polynomial = 1.0;
            break;
        case 1:
            hermite_polynomial = 2.0 * x;
            break;
        case 2:
            hermite_polynomial = 4.0 * pow(x, 2) - 2;
            break;
        case 3:
            hermite_polynomial = 8.0 * pow(x, 3) - 12 * x;
            break;
        case 4:
            hermite_polynomial = 16.0 * pow(x, 4) - 48.0 * pow(x, 2) + 12;
            break;
        case 5:
            hermite_polynomial = 32 * pow(x, 5) - 160 * pow(x, 3) + 120 * x;
            break;
        case 6:
            hermite_polynomial = 64 * pow(x, 6) - 480 * pow(x, 4) + 720 * pow(x, 2) - 120;
            break;
        case 7:
            hermite_polynomial = 128 * pow(x, 7) - 1344 * pow(x, 5) + 3360 * pow(x, 3) - 1680 * x;
            break;
        case 8:
            hermite_polynomial = 256 * pow(x, 8) - 3584 * pow(x, 6) + 13440 * pow(x, 4) - 13440 * pow(x, 2) + 1680;
            break;
        case 9:
            hermite_polynomial = 512 * pow(x, 9) - 9216 * pow(x, 7) + 48384 * pow(x, 5) - 80640 * pow(x, 3) + 30240 * x;
            break;
        case 10:
            hermite_polynomial = 1024 * pow(x, 10) - 23040 * pow(x, 8) + 161280 * pow(x, 6) - 403200 * pow(x, 4) + 302400 * pow(x, 2) - 30240;
            break;
        default:
            std::cout << "Error. Analytic expressions for the Hermite Polynomial of this degree has not been implemented.";
            exit(1);
            break;
    }
    return hermite_polynomial;
}
