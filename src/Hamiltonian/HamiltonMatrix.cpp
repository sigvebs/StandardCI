/* 
 * File:   HamiltonMatrix.cpp
 * Author: sigve
 * 
 * Created on October 22, 2012, 9:33 PM
 */

#include "HamiltonMatrix.h"


//------------------------------------------------------------------------------
HamiltonMatrix::HamiltonMatrix(Config *cfg, vector<bitset<BITS> > slaterDeterminants, mat interactions, vec spsEnergy) :
cfg(cfg), slaterDeterminants(slaterDeterminants), interactions(interactions), spsEnergy(spsEnergy)
{
    try {
        w = cfg->lookup("systemSettings.w");
        e0 = cfg->lookup("systemSettings.e0");
        wLaser = cfg->lookup("systemSettings.wLaser");
        wLaser *= w;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "SlaterDeterminants(Config *cfg, vector<vec> sps)::Error reading from 'systemSettings' object setting." << endl;
    }
#if DEBUG
    cout << "HamiltonMatrix(Config *cfg, vector<bitset<BITS> > slaterDeterminants, mat interactions, vec spsEnergy)" << endl;
    cout << "w = " << w << endl;
    cout << "e0 = " << e0 << endl;
    cout << "wLaser = " << wLaser << endl;
#endif
}
//------------------------------------------------------------------------------
void HamiltonMatrix::computeMatrixElements()
{
    cout << "Setting up the Hamilton matrix " << endl;
    int phase;
    int nStates = slaterDeterminants.size();

    bitset<BITS> newState;
    int i, j, k, l;
    double interactionElement;
    H = zeros(nStates, nStates);
    Ht = zeros(nStates, nStates);

    for (int st = 0; st < nStates; st++) {
        for (int b = 0; b < (int)interactions.n_rows; b++) {
            newState = slaterDeterminants[st];
            i = interactions(b, 0);
            j = interactions(b, 1);
            k = interactions(b, 2);
            l = interactions(b, 3);

            // Using the two body operator. Mathematics: a+(i)a+(j)a-(l)a-(k)|psi>
            phase = 1;

            newState = removeParticle(k, newState);
            phase *= sign(k, newState);

            newState = removeParticle(l, newState);
            phase *= sign(l, newState);

            newState = addParticle(j, newState);
            phase *= sign(j, newState);

            newState = addParticle(i, newState);
            phase *= sign(i, newState);

            if (newState[BITS - 1] != 1) {
                interactionElement = interactions(b, 4);
                // Searching for the new state to compute the matrix elements.
                for (int st2 = 0; st2 < (int)slaterDeterminants.size(); st2++) {
                    if (newState == slaterDeterminants[st2]) {
                        H(st, st2) += interactionElement * phase;
                        break;
                    }
                }
            }
        }

        // One body energies
        H(st, st) += spsEnergy(st);
    }

#if DEBUG
    cout << H << endl;

    //        // Symmetrizing the Hamilton matrix
    //        for (int i = 0; i < H.n_rows; i++) {
    //            for (int j = 0; j < i; j++) {
    //                H(j, i) = H(i, j);
    //            }
    //        }
    cout << "Completed generating the Hamilton matrix" << endl;
#endif
}
//------------------------------------------------------------------------------
void HamiltonMatrix::computeTimDepMatrixElements(int basisSize)
{
    cout << "Setting up the time dependent part of the Hamilton matrix." << endl;
    int phase;
    int nStates = slaterDeterminants.size();

    bitset<BITS> newState;
    double interactionElement = 1.0 / sqrt(2 * w);

    for (int st = 0; st < nStates; st++) {
        for (int p = 0; p < basisSize-1; p++) {

            // <Psi'|Psi^p_{p+1}>
            newState = slaterDeterminants[st];
            phase = 1;

            newState = removeParticle(p + 1, newState);
            phase *= sign(p + 1, newState);

            newState = addParticle(p, newState);
            phase *= sign(p, newState);

            if (newState[BITS - 1] != 1) {
                for (int st2 = 0; st2 < (int)slaterDeterminants.size(); st2++) {
                    if (newState == slaterDeterminants[st2]) {
                        Ht(st, st2) += sqrt(p+1)*interactionElement * phase;
                        break;
                    }
                }
            }

            // <Psi'|Psi^{p+1}_p>
            newState = slaterDeterminants[st];
            phase = 1;

            newState = removeParticle(p, newState);
            phase *= sign(p, newState);

            newState = addParticle(p + 1, newState);
            phase *= sign(p + 1, newState);

            if (newState[BITS - 1] != 1) {
                for (int st2 = 0; st2 < (int)slaterDeterminants.size(); st2++) {
                    if (newState == slaterDeterminants[st2]) {
                        Ht(st, st2) += sqrt(p+1)*interactionElement * phase;
                        break;
                    }
                }
            }
        }
    }

#if DEBUG
    cout << "computeTimDepMatrixElements(int basisSize)" << endl;
    cout << Ht << endl;
#endif
}
//------------------------------------------------------------------------------
int HamiltonMatrix::sign(int n, bitset<BITS> state) {
    int s = 1;
    for (int i = 0; i < n; i++) {
        if (state[i] != 0) {
            s *= -1;
        }
    }
    return s;
}
//------------------------------------------------------------------------------
bitset<BITS> HamiltonMatrix::addParticle(int n, bitset<BITS> state) {
    // Vacuum state
    if (state[BITS - 1])
        return state;

    bitset<BITS> a;
    a.set(n);
    bitset<BITS> comp = a & state;

    if (comp.count() == false) {
        state.set(n);
    } else {
        state.reset();
        state.set(BITS - 1);
    }

    return state;
}
//------------------------------------------------------------------------------
bitset<BITS> HamiltonMatrix::removeParticle(int n, bitset<BITS> state) {
    // Vacuum state
    if (state[BITS - 1])
        return state;

    bitset<BITS> a;
    a.set(n);
    bitset<BITS> comp = a & state;

    if (comp.count() == true) {
        state.set(n, 0);
    } else {
        state.reset();
        state.set(BITS - 1);
    }

    return state;
}
//------------------------------------------------------------------------------
mat HamiltonMatrix::getHamiltonian() {
    return H;
}
//------------------------------------------------------------------------------
mat HamiltonMatrix::operator()(double t) {
    return H + e0*sin(wLaser*t)*Ht;
}
//------------------------------------------------------------------------------
void HamiltonMatrix::setHamiltonian(mat h) {
    H = h;
    Ht = zeros(h.n_rows ,h.n_cols);
}
//------------------------------------------------------------------------------
mat HamiltonMatrix::evaluate(double t) {
    return H + e0*sin(wLaser*t)*Ht;
}
//------------------------------------------------------------------------------
int HamiltonMatrix::getDim() {
    return H.n_cols;
}
//------------------------------------------------------------------------------
