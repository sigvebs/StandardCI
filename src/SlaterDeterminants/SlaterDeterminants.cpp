/* 
 * File:   SlaterDeterminants.cpp
 * Author: sigve
 * 
 * Created on October 22, 2012, 7:53 PM
 */

#include "SlaterDeterminants.h"

//------------------------------------------------------------------------------

SlaterDeterminants::SlaterDeterminants() {
}

//------------------------------------------------------------------------------

/**
 * Constructor
 * @param nParticles number of particles
 * @param sps vector containing all orbitals
 */
SlaterDeterminants::SlaterDeterminants(Config *cfg, vector<vec> sps) : cfg(cfg), sps(sps) {
    
    try {
        nParticles = cfg->lookup("systemSettings.nParticles");
        conservedEigenSpin = cfg->lookup("systemSettings.conserveSpin");
        conservedEigenSpinValue = cfg->lookup("systemSettings.spinValue");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "SlaterDeterminants(Config *cfg, vector<vec> sps)::Error reading from 'systemSettings' object setting." << endl;
    }
#if DEBUG
    cout << "SlaterDeterminants(Config *cfg, vector<vec> sps)" << endl;
    cout << "nParticles = " << nParticles << endl;
    cout << "conservedEigenSpin = " << conservedEigenSpin << endl;
    cout << "conservedEigenSpinValue = " << conservedEigenSpinValue << endl;
#endif
}

//------------------------------------------------------------------------------

SlaterDeterminants::SlaterDeterminants(const SlaterDeterminants& orig) {
}

//------------------------------------------------------------------------------

SlaterDeterminants::~SlaterDeterminants() {
}


//------------------------------------------------------------------------------

/**
 * Creates all possible Slater determinants from a set of orbitals.
 */
void SlaterDeterminants::createSlaterDeterminants() {
    // Creating all possible states 
    int nSps = sps.size();
    vec state = zeros(nParticles, 1);

    // Creating an initial state
    for (int i = 0; i < nParticles; i++) {
        state(i) = i;
    }

    if (!conservedEigenSpin || checkEigenSpin(state)) {
        binStates.push_back(createBinaryState(state));
    }

    // Creating all possible slater determinants.
    while (true) {
        state = odometer(state, nSps, nParticles);

        if (sum(state) == 0)
            break;

        // Storing accepted states.
        if (!conservedEigenSpin || checkEigenSpin(state)) {
            binStates.push_back(createBinaryState(state));
        }
    }
}

//------------------------------------------------------------------------------

/**
 * Checks whether the eigenspin is conserved
 * @param state a slater determinant
 * @return Returns true if the eigenspin of a state is conserved.
 */
bool SlaterDeterminants::checkEigenSpin(vec state) {
    bool eigenSpin = true;
    int sumEigenSpin = 0;

    // Summing up the total eigenspin
    for (int i = 0; i < nParticles; i++) {
        sumEigenSpin += (sps[state(i)])[3];
    }

    // Checking if the total eigenspin of a state is conserved
    if (sumEigenSpin != conservedEigenSpinValue) {
        eigenSpin = false;
    }

    return eigenSpin;
}

//------------------------------------------------------------------------------

/**
 * Creates a binary representation of a slater determinant.
 * @param state a vector representing a slater determinant
 * @return returns a bit representation of a slater determinant
 */
bitset<BITS> SlaterDeterminants::createBinaryState(vec state) {
    // Constructing the binary representation
    bitset<BITS> binState;
    for (int i = 0; i < state.size(); i++) {
        try {
            binState.set(state(i));
        } catch (exception& e) {
            cout << "Exception: " << e.what() << endl;
            cout << "To run with this basis the number of bits must be increased. Current number of bits are " << BITS
                    << " , increase BITS > " << state(i) - 1 << ". Change BITS in defines.h and recompile." << endl;
            exit(1);
        }
    }

    return binState;
}

//------------------------------------------------------------------------------

/**
 * 
 * @param oldState
 * @param Nsp
 * @param N
 * @return 
 */
vec SlaterDeterminants::odometer(const vec &oldState, int Nsp, int N) {
    vec newState = oldState;
    double l;

    for (int j = N - 1; j >= 0; j--) {
        if (newState(j) < Nsp - N + j) {
            l = newState(j);
            for (int k = j; k < N; k++) {
                newState(k) = l + 1 + k - j;
            }
            return newState;
        }
    }

    newState = zeros(nParticles, 1);
    return newState;
}

//------------------------------------------------------------------------------

/**
 * 
 * @return returns a vector containing all Slater Determinants in binary form.
 */
vector<bitset<BITS> > SlaterDeterminants::getSlaterDeterminants() {
    return binStates;
}

//------------------------------------------------------------------------------

/**
 * 
 * @param spsEnergies
 * @return 
 */
vec SlaterDeterminants::computeSpsEnergies(vec spsEnergies) {
    int shells;
    double w;
    try {
        shells = cfg->lookup("systemSettings.shells");
        w = cfg->lookup("systemSettings.w");
 
    } catch (const SettingNotFoundException &nfex) {
        cerr << "SlaterDeterminants::computeSpsEnergies::Error reading from 'systemSettings' object setting." << endl;
    }
    
    int nStates = binStates.size();
    double E;
    for (int i = 0; i < nStates; i++) {
        E = 0;
        for (int j = 0; j < spsEnergies.n_elem; j++) {
            if (binStates[i][j]) {
                E += spsEnergies[j];
            }
        }
//        if (E > (shells + 1) * w)
//            binStates.erase(binStates.begin() + i);
    }

    // Calculating total sps energy
    vec energies = zeros(binStates.size(), 1);

    for (int i = 0; i < binStates.size(); i++) {
        for (int j = 0; j < spsEnergies.n_elem; j++) {
            if (binStates[i][j]) {
                energies[i] += spsEnergies[j];
            }
        }
    }
    return energies;
}

//------------------------------------------------------------------------------