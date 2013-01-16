/* 
 * File:   ForwardEuler.cpp
 * Author: sigve
 * 
 * Created on January 13, 2013, 10:15 PM
 */

#include "ForwardEuler.h"

//------------------------------------------------------------------------------

ForwardEuler::ForwardEuler() {
}

//------------------------------------------------------------------------------

ForwardEuler::ForwardEuler(Config *cfg, HamiltonMatrix *H, cx_vec C) : TimeIntegrator(cfg, H, C) {

}

//------------------------------------------------------------------------------

ForwardEuler::ForwardEuler(const ForwardEuler& orig) {
}

//------------------------------------------------------------------------------

void ForwardEuler::stepForward() {
    cx_vec CPrev = C;
    
    C.set_real(real(CPrev) + H->evaluate(t) * imag(CPrev) * dt);
    C.set_imag(imag(CPrev) - H->evaluate(t) * real(CPrev) * dt);
    
    t += dt;
}

//------------------------------------------------------------------------------

ForwardEuler::~ForwardEuler() {
}
//------------------------------------------------------------------------------
