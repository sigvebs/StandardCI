/* 
 * File:   MainApplication.cpp
 * Author: sigve
 *
 * Created on December 26, 2012, 3:26 PM
 */

#include "MainApplication.h"

//------------------------------------------------------------------------------
MainApplication::MainApplication(int* argc, char*** argv, string configFileName) : argc(argc), argv(argv) {
    cout << configFileName << endl;
    // Read the file. If there is an error, report it and exit.
    try {
        cfg.readFile(configFileName.c_str());
    } catch (const FileIOException &fioex) {
        cerr << "I/O error while reading config file." << endl;
        exit(EXIT_FAILURE);
    }
#if DEBUG
    cout << "Reading configuration from file" << endl;
#endif
}
//------------------------------------------------------------------------------
void MainApplication::runConfiguration() {

    bool cleanFiles = cfg.lookup("systemSettings.cleanFiles");
    if(cleanFiles)
        removeFiles();

    // Setting up a Harmonic Oscillator orbital basis
    mat interactionElements;
    vec orbitalElements;

    Basis basis = Basis(&cfg);
    basis.createBasis();

    string fileNameInteractionElements = createFileName("interactionElements");
    string fileNameOrbitalElements = createFileName("orbitalElements");

    // Setting the type of basis
    WaveFunction *wf;
    int dim = cfg.lookup("systemSettings.dim");

    switch(dim){
    case 1:
        wf = new HarmonicOscillator1d(&cfg);
        break;
    case 2:
        wf = new HarmonicOscillator2d(&cfg);
        break;
    }

    // Setting the integrator
    SpatialIntegrator *I;
    int spatialIntegrator;
    cfg.lookupValue("spatialIntegration.integrator", spatialIntegrator);

    switch (spatialIntegrator) {
    case MONTE_CARLO:
        I = new MonteCarloIntegrator(&cfg);
        break;
    case GAUSS_LAGUERRE:
        I = new GaussLaguerreIntegrator(&cfg, wf);
        break;
    case GAUSS_HERMITE:
        I = new GaussHermiteIntegrator(&cfg);
        break;
    case INTERACTION_INTEGRATOR:
        I = new InteractonIntegrator(&cfg);
        break;
    case MONTE_CARLO_IS:
        I = new MonteCarloImportanceSampled(&cfg);
        break;
    }

    cout << "spatialIntegrator = " << spatialIntegrator << endl;

    // Loading interaction elements and orbital elements from file
    // New elements are generated if needed.
    if (!interactionElements.load(fileNameInteractionElements)
            || !orbitalElements.load(fileNameOrbitalElements)) {
        cout << "Generating new orbital and interaction elements" << endl;
        basis.setIntegrator(I);
        basis.setWaveFunction(wf);
        basis.computeInteractionelements();
        basis.computeSpsEnergies();

        interactionElements = basis.getInteractionElements();
        orbitalElements = basis.getSpsEnergies();

        interactionElements.save(fileNameInteractionElements);
        orbitalElements.save(fileNameOrbitalElements);
    } else
        cout << "Interaction and orbital matrix elements found. Loading matrices from file." << endl;

    cout << interactionElements << endl;

    // Setting up all Slater Determinants
    vector<vec> states = basis.getStates();
    SlaterDeterminants SL = SlaterDeterminants(&cfg, states);
    SL.createSlaterDeterminants();
    vec spsEnergies = SL.computeSpsEnergies(orbitalElements);
    vector<bitset<BITS> > slaterDeterminants = SL.getSlaterDeterminants();
    cout << "Generated " << slaterDeterminants.size() << " Slater Determinants" << endl;

    // Setting up the Hamiltonian matrix
    HamiltonMatrix H = HamiltonMatrix(&cfg, slaterDeterminants, interactionElements, spsEnergies);
    mat h;
    string fileNameHamiltonian = createFileName("hamiltonian");

    if (!h.load(fileNameHamiltonian)) {
        cout << "Generating new Hamiltonian matrix" << endl;
        H.computeMatrixElements();
        h = H(0);
        h.save(fileNameHamiltonian);
    } else{
        cout << "Hamilton matrix found. Loading matrix from file." << endl;
        H.setHamiltonian(h);
    }

//    string fileNameImgHamiltonian = createFileName("imageHamiltonian");
//    mat hImage = 50*h;
//    hImage.save(fileNameImgHamiltonian, pgm_binary);
    

    // Finding the initial ground state
    vec eigval;
    mat eigvec;
    string fileNameEigval = createFileName("eigval");
    string fileNameEigvec = createFileName("eigvec");
    if (!eigval.load(fileNameEigval) || !eigvec.load(fileNameEigvec)) {
        cout << "Diagonalizing the Hamiltonian matrix" << endl;
        eig_sym(eigval, eigvec, h);
        eigval.save(fileNameEigval);
        eigvec.save(fileNameEigvec);
    } else
        cout << "Hamiltonian eigenvalues and eigenvectors loaded from file." << endl;

    // Printing results
    int samples = 0;

    switch (spatialIntegrator) {
    case MONTE_CARLO:
        cfg.lookupValue("spatialIntegration.MonteCarlo.samples", samples);
        break;
    case GAUSS_LAGUERRE:
        cfg.lookupValue("spatialIntegration.GaussLaguerre.samples", samples);
        break;
    case GAUSS_HERMITE:
        cfg.lookupValue("spatialIntegration.GaussHermite.samples", samples);
        break;
    }

    // Printing results
    double k = correlationFactor(eigvec.col(0));
    double w = cfg.lookup("systemSettings.w");
    cout << "\nIntegrator \t Shells \t Orbitals \t w \t Energy \t K \t Samples \n";
    cout << "-------------------------------------------------------------------------\n"
         << spatialIntegrator << " \t\t "
         << (int) cfg.lookup("systemSettings.shells") << " \t\t "
         << orbitalElements.size() << " \t\t "
         << w << " \t "
         << eigval.min() << " \t"
         << k << "\t"
         << samples << endl;
    cout << "-------------------------------------------------------------------------\n";


    // Time integration
    // Plotting resutls
    bool tIntegration;
    tIntegration = cfg.lookup("systemSettings.timeIntegration");
    if(tIntegration){
        string fileOverlap = createFileName("overlap");
        H.computeTimDepMatrixElements(orbitalElements.n_rows);
        int nCI = eigvec.n_rows;
        cx_vec C0(nCI);
        C0.set_real(eigvec.col(0));
        C0.set_imag(zeros(nCI, 1));
        TimeIntegrator *T;
        int tIntegrator = cfg.lookup("TimeIntegration.TimeIntegrator");

        switch (tIntegrator) {
        case FORWARD_EULER:
            T = new ForwardEuler(&cfg, &H, C0);
            break;
        case BACKWARD_EULER:
            cout << "BackwardEuler not yet implemented." << endl;
            break;
        case CRANK_NICOLSON:
            T = new CrankNicolson(&cfg, &H, C0);
            break;
        case EXPONENTIAL:
            T = new ExponentialPropagation(&cfg, &H, C0);
            break;
        }

        double wLaser = cfg.lookup("systemSettings.wLaser");
        wLaser *= w;
        double dt = cfg.lookup("TimeIntegration.dt");
        int end = floor(8*PI/(wLaser*dt));

        cx_vec C(nCI);
        cx_vec CPrev(nCI);
        vec overlap(end);
        C = C0;
        CPrev = C;

        for (int i = 0; i < end; i++) {
            T->stepForward();
            C = T->getCoefficients();
            overlap[i] = pow(abs(cdot(C, C0)), 2);
        }
        overlap.save(fileOverlap, arma_ascii);

        // Plotting resutls
        bool plotResult;
        plotResult = cfg.lookup("systemSettings.plotResult");
        if(plotResult)
            system("cd ..; cd PlottingResults; python plottingOverlap.py");
    }
    cout << "Run Complete" << endl;
}
//------------------------------------------------------------------------------
void MainApplication::finalize() {

}
//------------------------------------------------------------------------------
string MainApplication::createFileName(string baseName)
{
    string path;
    cfg.lookupValue("systemSettings.SIpath", path);
    return path + fName(baseName);
}
//------------------------------------------------------------------------------
string MainApplication::fName(string baseName)
{
    string fileName;
    //    string integrator;
    int integrator, shells, dim, nParticles, samples, basisType, coordinateType;
    double w, L;

    cfg.lookupValue("systemSettings.w", w);
    cfg.lookupValue("systemSettings.dim", dim);
    cfg.lookupValue("systemSettings.basisType", basisType);
    cfg.lookupValue("systemSettings.coordinateType", coordinateType);
    cfg.lookupValue("systemSettings.nParticles", nParticles);
    cfg.lookupValue("systemSettings.shells", shells);
    cfg.lookupValue("spatialIntegration.integrator", integrator);
    cfg.lookupValue("spatialIntegration.L", L);

    if (integrator == MONTE_CARLO)
        cfg.lookupValue("spatialIntegration.MonteCarlo.samples", samples);
    else if (integrator == GAUSS_LAGUERRE)
        cfg.lookupValue("spatialIntegration.GaussLaguerre.samples", samples);
    else if (integrator == GAUSS_HERMITE)
        cfg.lookupValue("spatialIntegration.GaussHermite.samples", samples);
    else if (integrator == INTERACTION_INTEGRATOR)
        samples = 0;
    else if (integrator == MONTE_CARLO_IS)
        cfg.lookupValue("spatialIntegration.MonteCarloIs.samples", samples);

    ostringstream convert;
    convert << "_basisType-" << basisType << "_coordinateType-"<< coordinateType << "_dim-" << dim << "_integrator-" << integrator << "_samples-" << samples << "_nParticles-" << nParticles << "_w-" << w << "_L-" << L << "_shells-" << shells ;

    fileName = baseName + convert.str() + ".mat";
    return fileName;
}
//------------------------------------------------------------------------------
double MainApplication::correlationFactor(vec p)
{
    double k = 0;
    for(int i=0;i<(int)p.n_elem; i++)
        k += pow(abs(p[i]),4);

    return 1.0/k;
}
//------------------------------------------------------------------------------
void MainApplication::removeFiles()
{
    string path;
    string command;

    cfg.lookupValue("systemSettings.SIpath", path);
    // command = "cd " + path + "; ls";
    command = "cd " + path + "; rm *" + fName("");
    cout << command << endl;
    system(command.c_str());
}
//------------------------------------------------------------------------------
