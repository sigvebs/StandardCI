/* 
 * File:   Basis.cpp
 * Author: sigve
 * 
 * Created on January 7, 2013, 4:22 PM
 */

#include "Basis.h"

//------------------------------------------------------------------------------ 

Basis::Basis() {
}

//------------------------------------------------------------------------------ 

Basis::Basis(Config *cfg) : cfg(cfg) {
    try {
        w = cfg->lookup("systemSettings.w");
        dim = cfg->lookup("systemSettings.dim");
        shells = cfg->lookup("systemSettings.shells");
        sIntegrator =  cfg->lookup("spatialIntegration.integrator");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Setting* systemSettings)::Error reading from 'systemSettings' object setting." << endl;
    }
    sqrtW = sqrt(w);
#if DEBUG
    cout << "Basis::Basis(Setting* systemSettings)" << endl;
    cout << "w = " << w << endl;
    cout << "dim = " << dim << endl;
#endif
}

//------------------------------------------------------------------------------ 

/**
 * Creates all possible orbitals in the specified number of shells.
 */
void Basis::createBasis() {
    // Generating all single particle states and energies        
    vec state(4);
    Orbital *orb;

    // TODO: generalize to different quantum numbers. 
    switch (dim) {
        case 1:
            for (int n = 0; n <= shells; n++) {
                for (int spin = 1; spin >= -1; spin -= 2) {
                    vec quantumNumbers(1);
                    quantumNumbers[0] = n;
                    //orb = new OrbitalHarmonicOscillator(systemSettings, quantumNumbers, spin);
                    orbitalStates.push_back(orb);
                    state(0) = states.size();
                    state(1) = n;
                    state(2) = 0;
                    state(3) = spin;
                    states.push_back(state);
                }
            }
            break;
        case 2:
            for (int k = 0; k <= shells; k++) {
                for (int i = 0; i <= floor(k / 2); i++) {
                    for (int l = -k; l <= k; l++) {
                        if ((2 * i + abs(l) + 1) == k) {

                            // Spin up
                            state(0) = states.size();
                            state(1) = i;
                            state(2) = l;
                            state(3) = 1;
                            states.push_back(state);

                            // Spin down
                            state(0) = states.size();
                            state(1) = i;
                            state(2) = l;
                            state(3) = -1;
                            states.push_back(state);
                        }
                    }
                }
            }
            break;
    }
    cout << states.size() << " orbitals created" << endl;
#if DEBUG
    cout << "void Basis::createBasis()" << endl;
    for (int i = 0; i < states.size(); i++) {
        for (int j = 0; j < states[0].n_elem; j++) {
            cout << states[i][j] << " ";
        }
        cout << endl;
    }
#endif
}
//------------------------------------------------------------------------------ 

/**
 * Computes all the interaction-elements between all possible configurations 
 * of orbitals.
 */
void Basis::computeInteractionelements() {
    cout << "Computing interaction elements" << endl;
    
    double tolerance = 1e-6;
    SpatialIntegrator *I;

    switch (sIntegrator) {
    case MONTE_CARLO:
       I = new MonteCarloIntegrator(cfg);
        break;
    case GAUSS_LAGUERRE:
         I = new GaussLaguerreIntegrator(cfg);
        break;
    case GAUSS_HERMITE:
         I = new GaussHermiteIntegrator(cfg);
        break;
    }

    int nStates = states.size();
    double E;

    vector<vec> interactionElements;
    vec interactionElement(5);

    for (int p = 0; p < nStates; p++) {
        for (int q = p + 1; q < nStates; q++) {
            for (int r = 0; r < nStates; r++) {
                for (int s = r + 1; s < nStates; s++) {
                    E = 0;

                    // Anti-symmetrized matrix elements
                    if (states[p][3] == states[r][3] && states[q][3] == states[s][3]) {
                        E += I->integrate(states[p][1], states[q][1], states[r][1], states[s][1]);
                    }

                    if (states[p][3] == states[s][3] && states[q][3] == states[r][3]) {
                        E -= I->integrate(states[p][1], states[q][1], states[s][1], states[r][1]);
                    }

                    if (abs(E) > tolerance) {
                        interactionElement << p << q << r << s << E;
                        interactionElements.push_back(interactionElement);
                    }
                }
            }
        }
    }
    cout << "Done " << endl;
    // Storing results in a matrix
    int nInteractionElements = interactionElements.size();
    intElements = zeros(nInteractionElements, 5);
    for (int i = 0; i < nInteractionElements; i++) {
        intElements.row(i) = trans(interactionElements[i]);
    }

#if DEBUG
    cout << "void Basis::computeInteractionelements()" << endl;
    cout << "InteractionELements = " << interactionElements.size() << endl;
#endif
}

//------------------------------------------------------------------------------ 

mat Basis::getInteractionElements() {
    return intElements;
}
//------------------------------------------------------------------------------ 

void Basis::computeSpsEnergies() {
    cout << "Computing orbital energies" << endl;
    double E;
    int nStates = states.size();
    spsEnergies = zeros(nStates, 1);

    vector<vec> spsElements;
    vec spsElement(3);

    switch (dim) {
        case 1:
            // Computing the one body operators
            for (int i = 0; i < states.size(); i++) {
                // Sps energies
                spsEnergies[i] = w * ((states[i])[1] + 0.5);
            }
            break;
    }

}

//------------------------------------------------------------------------------ 

vec Basis::getSpsEnergies() {
    return spsEnergies;
}
//------------------------------------------------------------------------------ 

vector<vec> Basis::getStates(){
    return states;
}

//------------------------------------------------------------------------------ 

Basis::~Basis() {
    // Cleaning up
    for (int i = 0; i < orbitalStates.size(); i++) {
        delete orbitalStates[i];
    }
}
