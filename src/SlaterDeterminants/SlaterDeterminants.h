/* 
 * File:   SlaterDeterminants.h
 * Author: sigve
 *
 * Created on October 22, 2012, 7:53 PM
 */

#ifndef SLATERDETERMINANTS_H
#define	SLATERDETERMINANTS_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <vector>
#include <map>
#include <bitset>
#include <libconfig.h++>
#include "../headers.h"

using namespace libconfig;
using namespace std;
using namespace arma;

class SlaterDeterminants {
public:
    SlaterDeterminants(Config *cfg, vector<vec> sps);
    void createSlaterDeterminants();
    vec odometer(const vec &, int, int);
    bool checkEigenSpin(vec state);
    bitset<BITS> createBinaryState(vec state);
    vector<bitset<BITS> > getSlaterDeterminants();
    vec computeSpsEnergies(vec spsEnergies);
private:
    Config *cfg;
    int nParticles;
    vector<vec> sps;
    vector<bitset<BITS> > binStates;
    bool conservedEigenSpin;
    int conservedEigenSpinValue;
};

#endif	/* SLATERDETERMINANTS_H */

