#include "testing.h"

//------------------------------------------------------------------------------
Testing::Testing()
{
    cout << "Configuring tests" << endl;

    // Reading configuration file - possibly change to TestConfig.cfg.
    string cfgFileName = "../testConfig.cfg";
    try {
        cfg.readFile(cfgFileName.c_str());
    } catch (const FileIOException &fioex) {
        cerr << "I/O error while reading config file." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
void Testing::runTests()
{
    cout << "Running tests" << endl;

    // Testing spatial integrators
    compareSpatialIntegrators();
}
//------------------------------------------------------------------------------
bool Testing::compareSpatialIntegrators()
{
    cout << "Testing integrators" << endl;

    SpatialIntegrator *MonteCarlo;
    SpatialIntegrator *GaussHermite;
    SpatialIntegrator *GaussLaguerre;
    SpatialIntegrator *MonteCarloIs;

    int dim = 1;
    double w = 0.25;

    cfg.lookup("systemSettings.dim") = dim;
    cfg.lookup("systemSettings.w") = w;


    WaveFunction *wf = new HarmonicOscillator1d(&cfg);

    // Test configurations
    int nQuantum = dim + 2;
    vec p = zeros(nQuantum), q = zeros(nQuantum), r = zeros(nQuantum), s = zeros(nQuantum);
    p << 0 << 1 << 0 << 0;
    q << 1 << -1 << 0 << 0;
    r << 2 << 1 << 0 << 0;
    s << 3 << -1 << 0 << 0;

    MonteCarlo      = new MonteCarloIntegrator(&cfg);
    GaussHermite    = new GaussHermiteIntegrator(&cfg);
    GaussLaguerre   = new GaussLaguerreIntegrator(&cfg, wf);
    MonteCarloIs    = new MonteCarloImportanceSampled(&cfg);

    double intMonteCarlo        = MonteCarlo->integrate(p, q, r, s);
    double intGaussHermite      = GaussHermite->integrate(p, q, r, s);
    double intGaussLaguerree    = GaussLaguerre->integrate(p, q, r, s);
    double intMonteCarloIs      = MonteCarloIs->integrate(p, q, r, s);

    // Printing results
    cout << "\n-------------------------------------------------------------\n";
    cout << "intMonteCarlo = " << intMonteCarlo << endl;
    cout << "intGaussHermite = " << intGaussHermite << endl;
    cout << "intGaussLaguerree = " << intGaussLaguerree << endl;
    cout << "intMonteCarloIs = " << intMonteCarloIs << endl;

    cout << endl;

    for(int i=1000; i<1e7; i *=10){
        cfg.lookup("spatialIntegration.MonteCarlo.samples") = i;
        MonteCarlo = new MonteCarloIntegrator(&cfg);
        intMonteCarlo = MonteCarlo->integrate(p, q, r, s);
        cout << i << "\t\t intMonteCarlo = " << intMonteCarlo << endl;
    }

    cout << endl;

    for(int i=1000; i<1e8; i *=10){
        cfg.lookup("spatialIntegration.MonteCarloIs.samples") = i;
        MonteCarloIs = new MonteCarloImportanceSampled(&cfg);
        intMonteCarloIs = MonteCarloIs->integrate(p, q, r, s);
        cout << i << "\t\t intMonteCarloIs = " << intMonteCarloIs << endl;
    }

    cout << endl;

    for(int i=80; i < 210; i +=20){
        cfg.lookup("spatialIntegration.GaussHermite.samples") = i;
        GaussHermite    = new GaussHermiteIntegrator(&cfg);
        intGaussHermite = GaussHermite->integrate(p, q, r, s);
        cout << i << "\t\t GaussHermite = " << intGaussHermite << endl;
    }

    cout << endl;

    for(int i=80; i < 200; i +=20){
        cfg.lookup("spatialIntegration.GaussLaguerre.samples") = i;
        GaussLaguerre    = new GaussLaguerreIntegrator(&cfg, wf);
        intGaussLaguerree = GaussLaguerre->integrate(p, q, r, s);
        cout << i << "\t\t intGaussLaguerree = " << intGaussLaguerree << endl;
    }
}
