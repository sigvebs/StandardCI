/* 
 * File:   MainApplication.h
 * Author: sigve
 *
 * Created on December 26, 2012, 3:26 PM
 */

#ifndef MAINAPPLICATION_H
#define	MAINAPPLICATION_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include <sstream>

#include "Basis/Basis.h"


#include "WaveFunction/wavefunction.h"
#include "WaveFunction/Implementations/harmonicoscillator1d.h"
#include "WaveFunction/Implementations/harmonicoscillator2d.h"

#include "Basis/SpatialIntegrator/SpatialIntegrator.h"
#include "Basis/SpatialIntegrator/MonteCarloIntegrator.h"
#include "Basis/SpatialIntegrator/GaussLaguerreIntegrator.h"
#include "Basis/SpatialIntegrator/GaussHermiteIntegrator.h"
#include "Basis/SpatialIntegrator/interactonintegrator.h"
#include "Basis/SpatialIntegrator/montecarloimportancesampled.h"

#include "SlaterDeterminants/SlaterDeterminants.h"
#include "Hamiltonian/HamiltonMatrix.h"

#include "TimeIntegrator/TimeIntegrator.h"
#include "TimeIntegrator/Integrator/ForwardEuler.h"
#include "TimeIntegrator/Integrator/CrankNicolson.h"
#include "TimeIntegrator/Integrator/exponentialpropagation.h"

using namespace std;
using namespace libconfig;

class MainApplication {
public:
    MainApplication(int* argc, char ***argv, string configFileName = "config.cfg");

    void runConfiguration();
    void finalize();
    string createFileName(string baseName);
    string fName(string baseName);
    double correlationFactor(vec p);
    void removeFiles();
    string results();
private:
    int* argc;
    char*** argv;
    Config cfg;
};

#endif	/* MAINAPPLICATION_H */

