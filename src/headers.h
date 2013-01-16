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

    //#define DEBUG 0
#define BITS 64
#define BUF_SIZE 33
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164

enum tIntegrator {
    FORWARD_EULER, BACKWARD_EULER, CRANK_NICOLSON
};
enum sIntegrator {
    MONTE_CARLO, GAUSS_LAGUERRE, GAUSS_HERMITE
};

#ifdef	__cplusplus
}
#endif

#endif	/* HEADERS_H */

