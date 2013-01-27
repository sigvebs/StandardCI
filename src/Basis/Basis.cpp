/* 
 * File:   Basis.cpp
 * Author: sigve
 *
 * Created on January 7, 2013, 4:22 PM
 */

#include "Basis.h"

//------------------------------------------------------------------------------ 
Basis::Basis(Config *cfg) : cfg(cfg)
{
    try {
        dim             = cfg->lookup("systemSettings.dim");
        coordinateType  = cfg->lookup("systemSettings.coordinateType");
        shells          = cfg->lookup("systemSettings.shells");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis::Basis(Setting* systemSettings)::Error reading from 'systemSettings' object setting." << endl;
    }

#if DEBUG
    cout << "Basis::Basis(Setting* systemSettings)" << endl;
    cout << "dim = " << dim << endl;
    cout << "coordinateType = " << coordinateType << endl;
    cout << "shells = " << shells << endl;
#endif
}
//------------------------------------------------------------------------------ 
void Basis::createBasis()
{
    switch(coordinateType){
    case CARTESIAN:
        createCartesianBasis();
        break;
    case POLAR:
        createPolarBasis();
        break;
    }
}
//------------------------------------------------------------------------------
void Basis::createCartesianBasis()
{
    // Generating all single particle states and energies
    vec state = zeros(dim + 2);

    switch (dim) {
    case 1:
        for (int n = 0; n <= shells; n++) {
            for (int spin = 1; spin >= -1; spin -= 2) {
                vec quantumNumbers(1);
                quantumNumbers[0] = n;
                state(0) = states.size();
                state(1) = spin;
                state(2) = n;
                states.push_back(state);
            }
        }
        break;
    case 2:
        for (int s = 0; s <= shells; s++) {
            for (int nx = 0; nx <= s; nx++) {
                for (int ny = 0; ny <= s; ny++) {
                    if(nx + ny == s){
                        for (int spin = 1; spin >= -1; spin -= 2) {
                            state(0) = states.size();
                            state(1) = spin;
                            state(2) = nx;
                            state(3) = ny;
                            states.push_back(state);
                        }
                    }
                }
            }
        }
        break;
    }
    cout << states.size() << " orbitals created" << endl;
#if DEBUG
    cout << "void Basis::createBasis()" << endl;
    for (int i = 0; i < (int)states.size(); i++) {
        for (int j = 0; j < (int)states[0].n_elem; j++) {
            cout << states[i][j] << " ";
        }
        cout << endl;
    }
#endif
}
//------------------------------------------------------------------------------
void Basis::createPolarBasis()
{
    vec state(dim +2);
    switch (dim) {
    case 1:
        cerr << "A polar basis in 1d does not make sense..." << endl;
        exit(1);
        break;
    case 2:
        for (int k = 0; k <= shells; k++) {
            for (int i = 0; i <= floor(k / 2); i++) {
                for (int l = -k; l <= k; l++) {
                    if ((2 * i + abs(l) + 1) == k) {

                        // Spin up
                        state(0) = states.size();
                        state(1) = 1;
                        state(2) = i;
                        state(3) = l;
                        states.push_back(state);

                        // Spin down
                        state(0) = states.size();
                        state(1) = -1;
                        state(2) = i;
                        state(3) = l;
                        states.push_back(state);
                    }
                }
            }
        }
        break;
    }
    cout << states.size() << " orbitals created" << endl;
#if DEBUG
    cout << "void Basis::createPolarBasis()" << endl;
    for (int i = 0; i < states.size(); i++) {
        for (int j = 0; j < states[0].n_elem; j++) {
            cout << states[i][j] << " ";
        }
        cout << endl;
    }
#endif
}
//------------------------------------------------------------------------------ 
void Basis::computeInteractionelements()
{
    cout << "Computing interaction elements" << endl;

    double tolerance = 1e-6;

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
                    if (states[p][1] == states[r][1] && states[q][1] == states[s][1]) {
                        E += I->integrate(states[p], states[q], states[r], states[s]);
                    }

                    if (states[p][1] == states[s][1] && states[q][1] == states[r][1]) {
                        E -= I->integrate(states[p], states[q], states[s], states[r]);
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
void Basis::setIntegrator(SpatialIntegrator *I)
{
    this->I = I;
}
//------------------------------------------------------------------------------
void Basis::setWaveFunction(WaveFunction *wf)
{
    this->wf = wf;
}
//------------------------------------------------------------------------------ 
mat Basis::getInteractionElements()
{
    return intElements;
}
//------------------------------------------------------------------------------ 
void Basis::computeSpsEnergies()
{
    cout << "Computing orbital energies" << endl;
    int nStates = states.size();
    spsEnergies = zeros(nStates, 1);

    // Computing the one body operators
    for (int i = 0; i < (int)states.size(); i++) {
        spsEnergies[i] = wf->getEnergy(states[i]);
    }
#if DEBUG
    cout << "void Basis::computeSpsEnergies()" << endl;
    cout << "spsEnergies = " << spsEnergies << endl;
#endif
}
//------------------------------------------------------------------------------ 
vec Basis::getSpsEnergies()
{
    return spsEnergies;
}
//------------------------------------------------------------------------------ 
vector<vec> Basis::getStates()
{
    return states;
}
//------------------------------------------------------------------------------ 
