/* 
 * File:   HamiltonMatrix.h
 * Author: sigve
 *
 * Created on October 22, 2012, 9:33 PM
 */

#ifndef HAMILTONMATRIX_H
#define	HAMILTONMATRIX_H

#include<iostream>
#include<fstream>
#include <vector>
#include <armadillo>

// For writing to file
#include <iostream>
#include <fstream>
#include <libconfig.h++>
#include "../headers.h"
#include <bitset>

using namespace libconfig;
using namespace std;
using namespace arma;

class HamiltonMatrix {
public:
    HamiltonMatrix();
    HamiltonMatrix(Config *cfg, vector<bitset<BITS> >, mat, vec);
    HamiltonMatrix(const HamiltonMatrix& orig);

    
    void computeMatrixElements();
    void computeTimDepMatrixElements(int basisSize);
    void setHamiltonian(mat h);
    
    // Binary operations
    bitset<BITS> addParticle(int, bitset<BITS>);
    bitset<BITS> removeParticle(int, bitset<BITS>);
    int sign(int, bitset<BITS> state);

    mat getHamiltonian();
    int getDim();
    mat evaluate(double t);
    mat operator()(double t);
    
    vec getSingleEnergy();
private:
    Config *cfg;
    mat H;
    mat Ht;
    vector<bitset<BITS> > slaterDeterminants;
    mat interactions;
    vec spsEnergy;
    double w;  // Frequency of the Harmonic oscillator.
    double e0; // Amplitude of laser.
    double wLaser; // Frequency of laser
};

#endif	/* HAMILTONMATRIX_H */

