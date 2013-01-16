/* 
 * File:   OrbitalHarmonicOscillator.h
 * Author: sigve
 *
 * Created on December 28, 2012, 3:32 PM
 */

#ifndef ORBITALHARMONICOSCILLATOR_H
#define	ORBITALHARMONICOSCILLATOR_H
#include "Orbital.h"

class OrbitalHarmonicOscillator: public Orbital {
public:
    OrbitalHarmonicOscillator();
    OrbitalHarmonicOscillator(Setting* systemSettings, vec qN, int sp);
    OrbitalHarmonicOscillator(const OrbitalHarmonicOscillator& orig);
    virtual ~OrbitalHarmonicOscillator();
    virtual double evaluate(vec x);
private:
    double hermitePolynomial(const int degree, const double x);
    vec quantumNumbers;
    int spin;
    int dim;
    double w;
    double sqrtW;
};

#endif	/* ORBITALHARMONICOSCILLATOR_H */

