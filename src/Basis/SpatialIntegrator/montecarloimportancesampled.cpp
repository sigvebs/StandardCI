#include "montecarloimportancesampled.h"

//------------------------------------------------------------------------------
MonteCarloImportanceSampled::MonteCarloImportanceSampled(Config *cfg) : SpatialIntegrator(cfg) {

    try {
        cfg->lookupValue("spatialIntegration.MonteCarloIs.samples", mcSamples);
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
double MonteCarloImportanceSampled::integrate(const vec &p, const vec &q, const
                                              vec &r, const vec &s) {
    double I; double x, theta, cosTheta, sinTheta;

    I = 0;
    for (int i = 0; i < mcSamples; i++) {
        x = sqrt(-log(1.0 - ran3(&idum)));
        theta = 2*PI*ran3(&idum);

        cosTheta = cos(theta);
        sinTheta = sin(theta);

        I +=      hermitePolynomial(p[2], x*cosTheta)
                * hermitePolynomial(q[2], x*sinTheta)
                * hermitePolynomial(r[2], x*cosTheta)
                * hermitePolynomial(s[2], x*sinTheta)
                / (sqrt( x*x*(1 - sin(2*theta))/w  + aa));
    }

    // The PI-factor has been reomved due to the change of variables.
    I *= 1 / (sqrt(pow(2, p[2] + q[2] + r[2] + s[2])
            * factorial(p[2])
            * factorial(q[2])
            * factorial(r[2])
            * factorial(s[2])));

    I /= mcSamples;

#if DEBUG
    cout << "p = " << p[2] << " q = " << q[2] << " r = " << r[2] << " s = " << s[2] << endl;
    cout << "MonteCarloIntegrator \t= " << I << "\t L = " << L << "\t MC samples = " << mcSamples << endl;
#endif
    return I;
}
//------------------------------------------------------------------------------
