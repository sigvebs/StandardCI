#ifndef TESTING_H
#define TESTING_H

#include <cstdlib>
#include <armadillo>
#include <libconfig.h++>

#include "src/Basis/SpatialIntegrator/SpatialIntegrator.h"
#include "src/Basis/SpatialIntegrator/GaussLaguerreIntegrator.h"
#include "src/Basis/SpatialIntegrator/GaussHermiteIntegrator.h"
#include "src/Basis/SpatialIntegrator/MonteCarloIntegrator.h"
#include "src/Basis/SpatialIntegrator/montecarloimportancesampled.h"

#include "src/WaveFunction/wavefunction.h"
#include "src/WaveFunction/Implementations/harmonicoscillator1d.h"

using namespace libconfig;
using namespace std;
using namespace arma;

class Testing
{
public:
    Testing();
    void runTests();
private:
    bool compareSpatialIntegrators();
    Config cfg;
};

#endif // TESTING_H
