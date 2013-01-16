/* 
 * File:   Orbital.cpp
 * Author: sigve
 * 
 * Created on December 28, 2012, 3:24 PM
 */
#include "Orbital.h"
//------------------------------------------------------------------------------

Orbital::Orbital() {
}
//------------------------------------------------------------------------------

Orbital::Orbital(vec quantumNumbers):quantumNumbers(quantumNumbers){
// TODO: does not seem to work this way 
#if DEBUG
    cout << "Orbital::Orbital(vec quantumNumbers):" << endl;
    cout << "quantumNumbers = " << quantumNumbers << endl;
#endif
}

//------------------------------------------------------------------------------

Orbital::Orbital(const Orbital& orig) {
}

//------------------------------------------------------------------------------

Orbital::~Orbital() {
}
//------------------------------------------------------------------------------

