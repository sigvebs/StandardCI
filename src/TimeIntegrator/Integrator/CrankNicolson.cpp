/* 
 * File:   CrankNicolson.cpp
 * Author: sigve
 * 
 * Created on January 14, 2013, 10:06 AM
 */

#include "CrankNicolson.h"

//------------------------------------------------------------------------------
CrankNicolson::CrankNicolson(Config *cfg, HamiltonMatrix *H, cx_vec C) : TimeIntegrator(cfg, H, C)
{
    n = H->getDim();
    I = eye<cx_mat>(n,n);
}
//------------------------------------------------------------------------------
void CrankNicolson::stepForward()
{
    cx_vec CPrev = C;
   
    H1 = cx_mat(eye(n,n), 0.5*dt*H->evaluate(t+dt));
    H2 = cx_mat(eye(n,n), -0.5*dt*H->evaluate(t)); 
    
    H1 = inv(H1);
    
    C = H1*H2*CPrev;

    t += dt;
}
//------------------------------------------------------------------------------
