/* 
 * File:   headers.h
 * Author: sigve
 *
 * Created on January 9, 2013, 6:41 PM
 */

#ifndef HEADERS_H
#define	HEADERS_H

#ifdef	__cplusplus
extern "C" {
#endif

#define DEBUG 0
#define TESTING 0
#define BITS 63
#define BUF_SIZE 33
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164

enum tIntegrator_t {
    FORWARD_EULER, BACKWARD_EULER, CRANK_NICOLSON, EXPONENTIAL
};
enum sIntegrator_t {
    MONTE_CARLO, GAUSS_LAGUERRE, GAUSS_HERMITE, INTERACTION_INTEGRATOR, MONTE_CARLO_IS
};
enum coordinateTypes_t {
    CARTESIAN, POLAR
};

#ifdef	__cplusplus
}
#endif

#endif	/* HEADERS_H */

