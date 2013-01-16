/* 
 * File:   Basis.h
 * Author: sigve
 *
 * Created on January 7, 2013, 4:22 PM
 */

#ifndef BASIS_H
#define	BASIS_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <vector>
#include <map>
#include <libconfig.h++>

#include "../headers.h"
#include "SpatialIntegrator/SpatialIntegrator.h"
#include "SpatialIntegrator/MonteCarloIntegrator.h"
#include "SpatialIntegrator/GaussLaguerreIntegrator.h"
#include "SpatialIntegrator/GaussHermiteIntegrator.h"
#include "../includes/lib.h"
#include "../Orbital/Orbital.h"
#include "../Orbital/OrbitalHarmonicOscillator.h"

using namespace libconfig;
using namespace std;
using namespace arma;

class Basis {
public:
    Basis();
    Basis(Config *cfg);
    Basis(const Basis& orig);
    virtual ~Basis();

    void createBasis();
    
    void computeSpsEnergies();
    vec getSpsEnergies();
    
    void computeInteractionelements();    
    mat getInteractionElements();
    
    vector<vec> getStates();
    
    // Wavefunction specifics
    double hermitePolynomial(const int degree, const double x);
    double waveFunction(double x, int n);
    double integrator(int p, int q, int r, int s);
    double expX(int p, int q);
    
private:
    //Setting* systemSettings;
    Config *cfg;

    string basisName;
    
    int sIntegrator;
    int dim;
    int shells;
    double w;
    double sqrtW;
    int maxRange;
    
    vector<vec> states;
    vector<Orbital*> orbitalStates;
    mat intElements;
    vec spsEnergies;
    
    int N;
};

#endif	/* BASIS_H */

