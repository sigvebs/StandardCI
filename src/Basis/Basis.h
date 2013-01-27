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
#include "../includes/lib.h"

#include "SpatialIntegrator/SpatialIntegrator.h"
#include "../WaveFunction/wavefunction.h"

using namespace libconfig;
using namespace std;
using namespace arma;

class Basis {
public:
    Basis(Config *cfg);

    void createBasis();
    void createCartesianBasis();
    void createPolarBasis();
    void computeSpsEnergies();
    void computeInteractionelements();
    void setIntegrator(SpatialIntegrator *I);
    void setWaveFunction(WaveFunction *wf);

    mat getInteractionElements();
    vec getSpsEnergies();
    vector<vec> getStates();
    
    // Wavefunction specifics
    double hermitePolynomial(const int degree, const double x);
    double waveFunction(double x, int n);
    double integrator(int p, int q, int r, int s);
    double expX(int p, int q);
    
protected:
    SpatialIntegrator *I;
    Config *cfg;

    string basisName;
    
    int sIntegrator;
    int coordinateType;

    int dim;
    int shells;

    int maxRange;
    
    vector<vec> states;
    WaveFunction *wf;
    mat intElements;
    vec spsEnergies;
    
    int N;
};

#endif	/* BASIS_H */

