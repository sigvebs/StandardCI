TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/MainApplication.cpp \
    src/Basis/Basis.cpp \
    src/Basis/SpatialIntegrator/SpatialIntegrator.cpp \
    src/Basis/SpatialIntegrator/MonteCarloIntegrator.cpp \
    src/Basis/SpatialIntegrator/GaussLaguerreIntegrator.cpp \
    src/Basis/SpatialIntegrator/GaussHermiteIntegrator.cpp \
    src/Hamiltonian/HamiltonMatrix.cpp \
    src/includes/ode.cpp \
    src/includes/lib.cpp \
    src/Integration/Integrator.cpp \
    src/Integration/GaussLagandre/gauss_legendre.c \
    src/Orbital/OrbitalHarmonicOscillator.cpp \
    src/Orbital/Orbital.cpp \
    src/SlaterDeterminants/SlaterDeterminants.cpp \
    src/TimeIntegrator/TimeIntegrator.cpp \
    src/TimeIntegrator/Integrator/ForwardEuler.cpp \
    src/TimeIntegrator/Integrator/CrankNicolson.cpp

OTHER_FILES += \
    src/includes/ode.sh \
    ../config.cfg

HEADERS += \
    src/MainApplication.h \
    src/headers.h \
    src/Basis/Basis.h \
    src/Basis/SpatialIntegrator/SpatialIntegrator.h \
    src/Basis/SpatialIntegrator/MonteCarloIntegrator.h \
    src/Basis/SpatialIntegrator/GaussLaguerreIntegrator.h \
    src/Basis/SpatialIntegrator/GaussHermiteIntegrator.h \
    src/Hamiltonian/HamiltonMatrix.h \
    src/includes/ode.hpp \
    src/includes/lib.h \
    src/Integration/Integrator.h \
    src/Integration/GaussLagandre/gauss_legendre.h \
    src/Orbital/OrbitalHarmonicOscillator.h \
    src/Orbital/Orbital.h \
    src/SlaterDeterminants/SlaterDeterminants.h \
    src/TimeIntegrator/TimeIntegrator.h \
    src/TimeIntegrator/Integrator/ForwardEuler.h \
    src/TimeIntegrator/Integrator/CrankNicolson.h


unix|win32: LIBS += -lconfig++ -larmadillo
