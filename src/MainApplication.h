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
#include "SlaterDeterminants/SlaterDeterminants.h"
#include "Hamiltonian/HamiltonMatrix.h"
#include "TimeIntegrator/TimeIntegrator.h"
#include "TimeIntegrator/Integrator/ForwardEuler.h"
#include "TimeIntegrator/Integrator/CrankNicolson.h"

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

