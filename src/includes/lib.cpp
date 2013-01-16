/*
 ** The library program  module
 **             lib.cpp
 **

TID time_step(int num)
 ** calculates start and stop time and returns the difference.                
 ** Input data:                            
 **    int number  = 1 for start - zero return values     
 **                = 2 for stop  - return time from last start
 **    TID run     - returns the time difference. 

void  **matrix(int row, int col, int num_bytes) 
 ** reserves dynamic memory for a two-dimensional matrix 
 ** using the C++ command new . No initialization of the elements. 
 ** Input data:                      
 **  int row      - number of  rows          
 **  int col      - number of columns        
 **  int num_bytes- number of bytes for each 
 **                 element                  
 ** Returns a void  **pointer to the reserved memory location.                                

void free_matrix(void **matr)
 ** releases the memory reserved by the function matrix() 
 ** for the two-dimensional matrix[][] 
 ** Input data:                          
 **  void  **matr - pointer to the matrix

void rk4(double *y, double *dydx, int n, double x, double h, double  *yout,
                    void (*derivs)(double, double *, double *))
 ** takes a set of variables y[1:n] for the function y(x) together with the
 ** reserves dynamic memory for a two-dimensional matrix 
 ** using the C++ command new . No initialization of the elements. 
 ** Input data:                      
 **  int row      - number of  rows          
 **  int col      - number of columns        
Æ*  int num_bytes- number of bytes for each 
 **                 element                  
 ** Returns a void  **pointer to the reserved memory location.                                

void free_matrix(void **matr)
 ** releases the memory reserved by the function matrix() for t    
 ** derivatives dydx[1:n] and uses the fourth-order Runge-Kutta method to
 ** advance the solution over an interval h and return incremented variables
 ** as yout[1:n], which not need to be a disstinct arra from y[1:n]. The
 ** users supply the routine derivs(x,y,dydx), which returns the derivatives
 ** dydx at x.

void ludcmp(double **a, int n, int *indx, double *d)
 ** takes as input a two-dimensional matrix a[][] of dimension n and
 ** replaces it by the LU decomposition of a rowwise permutation of
 ** itself. The results is stored in a[][] in the form given by 
 ** eq. (2.3.14) in "Numerical Recipe", sect. 2.3, page 45. The vector
 ** indx[] records the row permutation effected by the partial pivoting;
 ** d is output as +1 or -1 depending on whether the number of row
 ** interchanges was even or odd, respectively. This routine is used in
 ** combination with the function lubksb() to solve linear equations or
 ** invert a matrix. The function is slightly modified from the version
 ** in in Numerical recipe and uses memory allocation functions in the
 ** present module.

void lubksb(double **a, int n, int *indx, double *b)
 ** solves the set of linear equations A X = B of dimension n.
 ** a[][] is input, not as the matrix A[][] but rather as 
 ** its LU decomposition, determined by the function ludcmp(),
 ** indx[] is input as the permutation vector returned by 
 ** ludcmp(). b[] is input as the right-hand side vector B,
 ** The solution X is returned in B. The input data a[][],
 ** n and indx[] are not modified. This routine take into 
 ** account the possibility that b[] will begin with many
 ** zero elements, so it is efficient for use in matrix
 ** inversion.
 ** The function is slightly modified from the version in 
 ** in Numerical recipe.

void tqli(double d[], double e[], int n, double **z)
 ** determine eigenvalues and eigenvectors of a real symmetric
 ** tri-diagonal matrix, or a real, symmetric matrix previously
 ** reduced by function tred2[] to tri-diagonal form. On input,
 ** d[] contains the diagonal element and e[] the sub-diagonal
 ** of the tri-diagonal matrix. On output d[] contains the
 ** eigenvalues and  e[] is destroyed. If eigenvectors are
 ** desired z[][] on input contains the identity matrix. If
 ** eigenvectors of a matrix reduced by tred2() are required,
 ** then z[][] on input is the matrix output from tred2().
 ** On output, the k'th column returns the normalized eigenvector
 ** corresponding to d[k]. 
 ** The function is modified from the version in Numerical recipe.

void tred2(double **a, int n, double d[], double e[])
 ** perform a Housholder reduction of a real symmetric matrix
 ** a[][]. On output a[][] is replaced by the orthogonal matrix 
 ** effecting the transformation. d[] returns the diagonal elements
 ** of the tri-diagonal matrix, and e[] the off-diagonal elements, 
 ** with e[0] = 0.
 ** The function is modified from the version in Numerical recipe.

double pythag(double a, double b)
 ** The function is modified from the version in Numerical recipe.

void gauleg(double x1, double x2, double x[], double w[], int n)
 ** takes the lower and upper limits of integration x1, x2, calculates
 ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
 ** of length n of the Gauss--Legendre n--point quadrature formulae.

void jacobi(double** a, double* d, double** v, int n, int& nrot)
 ** Computes the eigenvalues and eigenvectors of the square symmetric matrix
 ** a  by use of the Jacobi method.
 ** It puts the eigenvalues in d and eigenvectors in v.
 ** n is an integer denoting the size of a
 **nrot keeps track of the number of rotations
 ** The function is as in the Numerical recipe

void jacobi_rot(double** a, double s, double tau, int i, int j, int k, int l)
 ** A helping function for jacobi making the actual rotations
 ** a is the matrix to be rotated, s is sine of the rotation
 ** angle, tau is s/(1 + c) where c is cosine of the angle.
 ** The integers i-l denotes the matrix element to be
 ** rotated.


double rectangle_rule(double a, double b, int n, double (*func)(double))
 ** integration by rectangle rule, in a and b and number of points n and name
 ** of function
double trapezoidal_rule(double a, double b, int n, double (*func)(double))
 ** integration by trapezoidal rule, in a and b and number of points n and name
 ** of function
void spline(double x[], double y[], int n, double yp1, double yp2, double y2[])
 ** takes as input x[0,..,n - 1] and y[0,..,n - 1] containing a tabulation
 ** y_i = f(x_i) with x_0 < x_1 < .. < x_(n - 1) together with yp_1 and yp2
 ** for first derivatives  f(x) at x_0 and x_(n-1), respectively. Then the
 ** function returns y2[0,..,n-1] which contanin the second derivatives of
 ** f(x_i)at each point x_i. If yp1 and/or yp2 is larger than the constant
 ** INFINITY the function will put corresponding second derivatives to zero.

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
 ** takes xa[0,..,n - 1] and y[0,..,n - 1] which tabulates a function 
 ** (with the xa[i]'s in order) and given ya[0,..,n - 1], which is the
 ** output from function spline() and with given value of x returns a 
 ** cubic--spline interpolation value y.

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
 ** takes as input xa[0,..,n-1] and ya[0,..,n-1] together with a given value
 ** of x and returns a value y and an error estimate dy. If P(x) is a polynomial
 ** of degree N - 1 such that P(xa_i) = ya_i, i = 0,..,n-1, then the returned 
 ** value is y = P(x). 

double rtbis(double (*func)(double), double x1, double x2, double xacc)
 ** calculates a root between x1 and x2 of a function
 ** pointed to by (*func) using the method of bisection  
 ** The root is returned with an accuracy of +- xacc.

double rtsec(double (*func)(double), double x1, double x2, double xacc)
 ** calculates a root between x1 and x2 of a function
 ** pointed to by (*func) using the secant method.
 ** The root is returned with an accuracy of +- xacc.

double rtnewt(void (*funcd)(double, double *, double *), double x1, double x2, double xacc)
 ** calculates a root between x1 and x2 of a function pointed to
 ** by (*funcd) using the Newton-Raphson method. The user-defined
 ** function funcd() returns both the function value and its first
 ** derivative at the point x,
 ** The root is returned with an accuracy of +- xacc.

double zbrent(double (*func)(double), double x1, double x2, double xacc)
 ** calculates a root between x1 and x2 of a function
 ** pointed to by (*funcd) using the Brent's method.
 ** The root is returned with an accuracy of +- xacc.

double ran0(long *idum)
 ** is an "Minimal" random number generator of Park and Miller
 ** (see Numerical recipe page 279). Set or reset the input value
 ** idum to any integer value (except the unlikely value MASK)
 ** to initialize the sequence; idum must not be altered between
 ** calls for sucessive deviates in a sequence.
 ** The function returns a uniform deviate between 0.0 and 1.0.

double ran1(long *idum)
 ** is an "Minimal" random number generator of Park and Miller
 ** (see Numerical recipe page 280) with Bays-Durham shuffle and
 ** added safeguards. Call with idum a negative integer to initialize;
 ** thereafter, do not alter idum between sucessive deviates in a
 ** sequence. RNMX should approximate the largest floating point value
 ** that is less than 1.
 ** The function returns a uniform deviate between 0.0 and 1.0
 ** (exclusive of end-point values).

double ran2(long *idum)
 ** is a long periode (> 2 x 10^18) random number generator of 
 ** L'Ecuyer and Bays-Durham shuffle and added safeguards.
 ** Call with idum a negative integer to initialize; thereafter,
 ** do not alter idum between sucessive deviates in a
 ** sequence. RNMX should approximate the largest floating point value
 ** that is less than 1.
 ** The function returns a uniform deviate between 0.0 and 1.0
 ** (exclusive of end-point values).
   
double ran3(long *idum)
 ** returns a uniform random number deviate between 0.0 and 1.0. Set
 ** the idum to any negative value to initialize or reinitialize the
 ** sequence. Any large MBIG, and any small (but still large) MSEED
 ** can be substituted for the present values. 
 */

#include "lib.h"

/*
 ** The function                           
 **      TID time_step(..)                 
 ** calculates start and stop time and returns the difference.                
 ** Input data:                            
 **    int number  = 1 for start - zero return values     
 **                = 2 for stop  - return time from last start
 **    TID run     - returns the time difference. 
 */

TID time_step(int num) {
    unsigned long long int
    num_sec;

    static long
    zsec = 0, zusec = 0;
    double
    retval;
    TID
    ex_time;
    struct timeval
    tp;

    if (num == 1) { // initialization of time

        zsec = tp.tv_sec;
        zusec = tp.tv_usec;

        ex_time.sec = 0;
        ex_time.min = 0;
        ex_time.hour = 0;
    } else if (num == 2) {

        retval = (double) (tp.tv_sec - zsec) + (tp.tv_usec - zusec) * 0.000001;

        num_sec = (unsigned long long int) retval;
        ex_time.sec = num_sec % 60;
        ex_time.min = num_sec / 60;
        ex_time.hour = ex_time.min / 60;
        ex_time.min = ex_time.min % 60;
    } else {
        printf("\n\nError in function time_step(): ");
        printf("\nInput data num = %d is wrong !!\n\n", num);
        exit(1);
    }
    return ex_time;

} // End: function time_step()

/*
 * The function                             
 *      void  **matrix()                    
 * reserves dynamic memory for a two-dimensional matrix 
 * using the C++ command new . No initialization of the elements. 
 * Input data:                      
 *  int row      - number of  rows          
 *  int col      - number of columns        
 *  int num_bytes- number of bytes for each 
 *                 element                  
 * Returns a void  **pointer to the reserved memory location.                                
 */

void **matrix(int row, int col, int num_bytes) {
    int i, num;
    char **pointer, *ptr;

    pointer = new(nothrow) char* [row];
    if (!pointer) {
        cout << "Exception handling: Memory allocation failed";
        cout << " for " << row << "row addresses !" << endl;
        return NULL;
    }
    i = (row * col * num_bytes) / sizeof (char);
    pointer[0] = new(nothrow) char [i];
    if (!pointer[0]) {
        cout << "Exception handling: Memory allocation failed";
        cout << " for address to " << i << " characters !" << endl;
        return NULL;
    }
    ptr = pointer[0];
    num = col * num_bytes;
    for (i = 0; i < row; i++, ptr += num) {
        pointer[i] = ptr;
    }

    return (void **) pointer;

} // end: function void **matrix()

/*
 * The function                         
 *      void free_matrix()              
 * releases the memory reserved by the function matrix() 
 *for the two-dimensional matrix[][] 
 * Input data:                          
 *  void far **matr - pointer to the matrix
 */

void free_matrix(void **matr) {

    delete [] (char *) matr[0];
    delete [] matr;

} // End:  function free_matrix() 

/*
 ** The function
 **              rk4()
 ** takes a set of variables y[1:n] for the function y(x) together with the
 ** derivatives dydx[1:n] and uses the fourth-order Runge-Kutta method to
 ** advance the solution over an interval h and return incremented variables
 ** as yout[1:n], which not need to be a disstinct arra from y[1:n]. The
 ** users supply the routine derivs(x,y,dydx), which returns the derivatives
 ** dydx at x.
 */

void rk4(double *y, double *dydx, int n, double x, double h, double *yout,
        void (*derivs)(double, double *, double *)) {
    int i;
    double xh, hh, h6, *dym, *dyt, *yt;


    // local memory allocation

    dym = new(nothrow) double [n];
    if (!dym) {
        printf("\n\nError in function rk4():");
        printf("\nNot enough memory for dym[%d]\n", n);
        exit(1);
    }
    dyt = new(nothrow) double [n];
    if (!dyt) {
        printf("\n\nError in function rk4():");
        printf("\nNot enough memory for dyt[%d]\n", n);
        exit(1);
    }
    yt = new(nothrow) double [n];
    if (!yt) {
        printf("\n\nError in function rk4():");
        printf("\nNot enough memory for yt[%d]\n", n);
        exit(1);
    }

    hh = h * 0.5;
    h6 = h / 6.0;
    xh = x + hh;

    for (i = 0; i < n; i++) { // first step
        yt[i] = y[i] + hh * dydx[i];
    }

    (*derivs)(xh, yt, dyt); // second step

    for (i = 1; i < n; i++) {
        yt[i] = y[i] + hh * dyt[i];
    }

    (*derivs)(xh, yt, dym); // third step 

    for (i = 1; i < n; i++) {
        yt[i] = y[i] + h * dym[i];
        dym[i] += dyt[i];
    }

    (*derivs)(x + h, yt, dyt); // fourth step 

    // acummulate increments with proper weights

    for (i = 0; i < n; i++) {
        yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    }

    delete [] yt; // release local memory
    delete [] dyt;
    delete [] dym;

} // End: function rk4()

/*
 ** The function
 **       ludcmp()
 ** takes as input a two-dimensional matrix a[][] of dimension n and
 ** replaces it by the LU decomposition of a rowwise permutation of
 ** itself. The results is stored in a[][] in the form given by 
 ** eq. (2.3.14) in "Numerical Recipe", sect. 2.3, page 45. The vector
 ** indx[] records the row permutation effected by the partial pivoting;
 ** d is output as +1 or -1 depending on whether the number of row
 ** interchanges was even or odd, respectively. This routine is used in
 ** combination with the function lubksb() to solve linear equations or
 ** invert a matrix. The function is slightly modified from the version
 ** in in Numerical recipe and uses memory allocation functions in the
 ** present module.
 */

void ludcmp(double **a, int n, int *indx, double *d) {
    int i, imax, j, k;
    double big, dum, sum, temp, *vv;

    vv = new(nothrow) double [n];
    if (!vv) {
        printf("\n\nError in function ludcm():");
        printf("\nNot enough memory for vv[%d]\n", n);
        exit(1);
    }

    *d = 1.0; // no row interchange yet
    for (i = 0; i < n; i++) { // loop over rows to get scaling information
        big = ZERO;
        for (j = 0; j < n; j++) {
            if ((temp = fabs(a[i][j])) > big) big = temp;
        }
        if (big == ZERO) {
            printf("\n\nSingular matrix in routine ludcmp()\n");
            exit(1);
        }
        vv[i] = 1.0 / big; // save scaling */
    } // end i-loop */

    for (j = 0; j < n; j++) { // loop over columns of Crout's method
        for (i = 0; i < j; i++) { // not i = j
            sum = a[i][j];
            for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = ZERO; // initialization for search for largest pivot element
        for (i = j; i < n; i++) {
            sum = a[i][j];
            for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        } // end i-loop
        if (j != imax) { // do we need to interchange rows ?
            for (k = 0; k < n; k++) { // yes
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            (*d) *= -1; // and change the parit of d
            vv[imax] = vv[j]; // also interchange scaling factor 
        }
        indx[j] = imax;
        if (fabs(a[j][j]) < ZERO) a[j][j] = ZERO;

        /*
         ** if the pivot element is zero the matrix is singular
         ** (at least to the precision of the algorithm). For 
         ** some application of singular matrices, it is desirable
         ** to substitute ZERO for zero,
         */

        if (j < (n - 1)) { // divide by pivot element 
            dum = 1.0 / a[j][j];
            for (i = j + 1; i < n; i++) a[i][j] *= dum;
        }
    } // end j-loop over columns

    delete [] vv; // release local memory

} // End: function ludcmp()

/*
 ** The function 
 **             lubksb()
 ** solves the set of linear equations A X = B of dimension n.
 ** a[][] is input, not as the matrix A[][] but rather as 
 ** its LU decomposition, determined by the function ludcmp(),
 ** indx[] is input as the permutation vector returned by 
 ** ludcmp(). b[] is input as the right-hand side vector B,
 ** The solution X is returned in B. The input data a[][],
 ** n and indx[] are not modified. This routine take into 
 ** account the possibility that b[] will begin with many
 ** zero elements, so it is efficient for use in matrix
 ** inversion.
 ** The function is slightly modified from the version in 
 ** in Numerical recipe.
 */

void lubksb(double **a, int n, int *indx, double *b) {
    int i, ii = -1, ip, j;
    double sum;

    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii > -1) for (j = ii; j < i; j++) sum -= a[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
} // End: function lubksb()

/*
 ** The function
 **                 tqli()
 ** determine eigenvalues and eigenvectors of a real symmetric
 ** tri-diagonal matrix, or a real, symmetric matrix previously
 ** reduced by function tred2[] to tri-diagonal form. On input,
 ** d[] contains the diagonal element and e[] the sub-diagonal
 ** of the tri-diagonal matrix. On output d[] contains the
 ** eigenvalues and  e[] is destroyed. If eigenvectors are
 ** desired z[][] on input contains the identity matrix. If
 ** eigenvectors of a matrix reduced by tred2() are required,
 ** then z[][] on input is the matrix output from tred2().
 ** On output, the k'th column returns the normalized eigenvector
 ** corresponding to d[k]. 
 ** The function is modified from the version in Numerical recipe.
 */

void tqli(double *d, double *e, int n, double **z) {
    register int m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;

    for (i = 1; i < n; i++) e[i - 1] = e[i];
    e[n] = 0.0;
    for (l = 0; l < n; l++) {
        iter = 0;
        do {
            for (m = l; m < n - 1; m++) {
                dd = fabs(d[m]) + fabs(d[m + 1]);
                if ((double) (fabs(e[m]) + dd) == dd) break;
            }
            if (m != l) {
                if (iter++ == 30) {
                    printf("\n\nToo many iterations in tqli.\n");
                    exit(1);
                }
                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p = 0.0;
                for (i = m - 1; i >= l; i--) {
                    f = s * e[i];
                    b = c * e[i];
                    e[i + 1] = (r = pythag(f, g));
                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;
                    for (k = 0; k < n; k++) {
                        f = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    } /* end k-loop */
                } /* end i-loop */
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            } /* end if-loop for m != 1 */
        } while (m != l);
    } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */

/*
 ** The function
 **                tred2()
 ** perform a Housholder reduction of a real symmetric matrix
 ** a[][]. On output a[][] is replaced by the orthogonal matrix 
 ** effecting the transformation. d[] returns the diagonal elements
 ** of the tri-diagonal matrix, and e[] the off-diagonal elements, 
 ** with e[0] = 0.
 ** The function is modified from the version in Numerical recipe.
 */

void tred2(double **a, int n, double *d, double *e) {
    register int l, k, j, i;
    double scale, hh, h, g, f;

    for (i = n - 1; i > 0; i--) {
        l = i - 1;
        h = scale = 0.0;
        if (l > 0) {
            for (k = 0; k <= l; k++)
                scale += fabs(a[i][k]);
            if (scale == 0.0) // skip transformation
                e[i] = a[i][l];
            else {
                for (k = 0; k <= l; k++) {
                    a[i][k] /= scale; // used scaled a's for transformation
                    h += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i] = scale*g;
                h -= f * g;
                a[i][l] = f - g;
                f = 0.0;

                for (j = 0; j <= l; j++) {
                    a[j][i] = a[i][j] / h; // can be omitted if eigenvector not wanted
                    g = 0.0;
                    for (k = 0; k <= j; k++) {
                        g += a[j][k] * a[i][k];
                    }
                    for (k = j + 1; k <= l; k++)
                        g += a[k][j] * a[i][k];
                    e[j] = g / h;
                    f += e[j] * a[i][j];
                }
                hh = f / (h + h);
                for (j = 0; j <= l; j++) {
                    f = a[i][j];
                    e[j] = g = e[j] - hh*f;
                    for (k = 0; k <= j; k++)
                        a[j][k] -= (f * e[k] + g * a[i][k]);
                }
            }
        }// end if-loop for l > 1
        else {
            e[i] = a[i][l];
        }
        d[i] = h;
    } // end i-loop
    d[0] = 0.0;
    e[0] = 0.0;

    /* Contents of this loop can be omitted if eigenvectors not
     ** wanted except for statement d[i]=a[i][i];
     */

    for (i = 0; i < n; i++) {
        l = i - 1;
        if (d[i]) {
            for (j = 0; j <= l; j++) {
                g = 0.0;
                for (k = 0; k <= l; k++) {
                    g += a[i][k] * a[k][j];
                }
                for (k = 0; k <= l; k++) {
                    a[k][j] -= g * a[k][i];
                }
            }
        }
        d[i] = a[i][i];
        a[i][i] = 1.0;
        for (j = 0; j <= l; j++) {
            a[j][i] = a[i][j] = 0.0;
        }
    }
} // End: function tred2(), (C) Copr. 1986-92 Numerical Recipes Software )

double pythag(double a, double b) {
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
    else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}
// End: function pythag(), (C) Copr. 1986-92 Numerical Recipes Software )%.

/*
 ** The function 
 **              gauleg()
 ** takes the lower and upper limits of integration x1, x2, calculates
 ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
 ** of length n of the Gauss--Legendre n--point quadrature formulae.
 */

void gauleg(double x1, double x2, double x[], double w[], int n) {
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;
    double const pi = 3.14159265359;
    double *x_low, *x_high, *w_low, *w_high;

    m = (n + 1) / 2; // roots are symmetric in the interval
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    x_low = x; // pointer initialization
    x_high = x + n - 1;
    w_low = w;
    w_high = w + n - 1;

    for (i = 1; i <= m; i++) { // loops over desired roots
        z = cos(pi * (i - 0.25) / (n + 0.5));

        /*
         ** Starting with the above approximation to the ith root
         ** we enter the mani loop of refinement bt Newtons method.
         */

        do {
            p1 = 1.0;
            p2 = 0.0;

            /*
             ** loop up recurrence relation to get the
             ** Legendre polynomial evaluated at x
             */

            for (j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }

            /*
             ** p1 is now the desired Legrendre polynomial. Next compute
             ** ppp its derivative by standard relation involving also p2,
             ** polynomial of one lower order.
             */

            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp; // Newton's method
        } while (fabs(z - z1) > ZERO);

        /* 
         ** Scale the root to the desired interval and put in its symmetric
         ** counterpart. Compute the weight and its symmetric counterpart
         */

        *(x_low++) = xm - xl * z;
        *(x_high--) = xm + xl * z;
        *w_low = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        *(w_high--) = *(w_low++);
    }
} // End_ function gauleg()

/*
 ** The function 
 **           jacobi_rot()
 ** A helping function for jacobi making the actual rotations
 ** a is the matrix to be rotated, s is sine of the rotation
 ** angle, tau is s/(1 + c) where c is cosine of the angle.
 ** The integers i-l denotes the matrix element to be
 ** rotated.
 **
 */

inline void jacobi_rot(double** a, double s, double tau, int i, int j, int k, int l) {
    double g, h;

    g = a[i][j];
    h = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);

}//End function jacobi_rot

/*
 ** The function 
 **              jacobi()
 ** Computes the eigenvalues and eigenvectors of the square symmetric matrix
 ** A  by use of the Jacobi method.
 ** It puts the eigenvalues in d and eigenvectors in v.
 ** n is an integer denoting the size of A
 ** The function is as in the Numerical recipe
 */
void jacobi(double** a, double* d, double** v, int n, int &nrot) {
    int i, j, ip, iq;
    double tresh, theta, tau, t, sm, s, h, g, c;

    double* b = new double[n];
    double* z = new double[n];
    for (ip = 0; ip < n; ip++) {
        for (iq = 0; iq < n; iq++) {
            v[ip][iq] = 0.0; //Initializing v to the identity matrix
            v[ip][ip] = 1.0;
        }
    }

    for (ip = 0; ip < n; ip++) {
        b[ip] = d[ip] = a[ip][ip]; //Initializing d and b to the diagonal of a
        z[ip] = 0.0; //z will accumulate terms of the form
        //t*a[ip][iq]
    }


    nrot = 0;
    for (i = 1; i <= 50; i++) {
        sm = 0.0;
        for (ip = 0; ip < n - 1; ip++) {
            for (iq = ip + 1; iq < n; iq++) {
                sm += fabs(a[ip][iq]); //Sum magnitude of off-diagonal elements
            }
        }
        if (sm == 0.0) {
            return; //The normal return at convergence
        }
        if (i < 4) {
            tresh = 0.2 * sm / (n * n); //On the first four sweeps
        } else {
            tresh = 0.0; //... thereafter
        }
        for (ip = 0; ip < n - 1; ip++) {
            for (iq = ip + 1; iq < n; iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                //After four sweeps we skip the rotation if the off-diagonal element is small
                if (i > 4 && (fabs(d[ip]) + g) == fabs(d[ip])
                        && (fabs(d[iq]) + g) == fabs(d[iq])) {
                    a[ip][iq] = 0.0;
                } else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if ((fabs(h) + g) == fabs(h)) {
                        t = (a[ip][iq]) / h;
                    } else {
                        theta = 0.5 * h / (a[ip][iq]);
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) {
                            t = -t;
                        }
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t*c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j < ip; j++) {
                        jacobi_rot(a, s, tau, j, ip, j, iq); // Rotations for 0 <= j < ip
                    }
                    for (j = ip + 1; j < iq; j++) {
                        jacobi_rot(a, s, tau, ip, j, j, iq); // Rotations for ip < j < iq
                    }
                    for (j = iq + 1; j < n; j++) {
                        jacobi_rot(a, s, tau, ip, j, iq, j); // Rotations for q < j < n
                    }
                    for (j = 0; j < n; j++) {
                        jacobi_rot(v, s, tau, j, ip, j, iq); //Updating v
                    }
                    nrot++;
                }
            }
        }
        for (ip = 0; ip < n; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    printf("\n\nToo many iterations in routine jacobi.\n");
    exit(1);
}//End function jacobi()

/*
 ** The function 
 **           spline()
 ** takes as input x[0,..,n - 1] and y[0,..,n - 1] containing a tabulation
 ** y_i = f(x_i) with x_0 < x_1 < .. < x_(n - 1) together with yp_1 and yp2
 ** for first derivatives  f(x) at x_0 and x_(n-1), respectively. Then the
 ** function returns y2[0,..,n-1] which contanin the second derivatives of
 ** f(x_i)at each point x_i. If yp1 and/or yp2 is larger than the constant
 ** INFINITY the function will put corresponding second derivatives to zero.
 */

void spline(double x[], double y[], int n, double yp1, double yp2, double y2[]) {
    int i, k;
    double p, qn, sig, un, *u;

    u = new(nothrow) double [n];
    if (!u) {
        printf("\n\nError in function spline():");
        printf("\nNot enough memory for u[%d]\n", n);
        exit(1);
    }

    if (yp1 > INFINITY) y2[0] = u[0] = 0.0;
    else {
        y2[0] = -0.5;
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    for (i = 1; i < (n - 1); i++) {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (yp2 > INFINITY) qn = un = ZERO;
    else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2])) * (yp2 - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

    for (k = n - 2; k >= 0; k--) {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }
    delete [] u;
    ;

} // End: function spline()

/*
 ** The function 
 **           splint()
 ** takes xa[0,..,n - 1] and y[0,..,n - 1] which tabulates a function 
 ** (with the xa[i]'s in order) and given ya[0,..,n - 1], which is the
 ** output from function spline() and with given value of x returns a 
 ** cubic--spline interpolation value y.
 */

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y) {
    int klo, khi, k;
    double h, b, a;

    klo = 0;
    khi = n - 1;
    while ((khi - klo) > 1) { // binary search
        k = (khi + klo) >> 1;
        if (xa[k] > x) khi = k;
        else klo = k;
    }
    h = xa[khi] - xa[klo];
    if (fabs(h) < ZERO) {
        printf("\n\n Error in function splint(): ");
        printf("\n The difference h = %4.1E -- too small\n", h);
        exit(1);
    }
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo]
            + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;

} // End: function splint()

/*
 ** The function
 **            polint()
 ** takes as input xa[0,..,n-1] and ya[0,..,n-1] together with a given value
 ** of x and returns a value y and an error estimate dy. If P(x) is a polynomial
 ** of degree N - 1 such that P(xa_i) = ya_i, i = 0,..,n-1, then the returned 
 ** value is y = P(x). 
 */

void polint(double xa[], double ya[], int n, double x, double *y, double *dy) {
    int i, m, ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;

    dif = fabs(x - xa[0]);

    c = new(nothrow) double [n];
    if (!c) {
        printf("\n\nError in function polint():");
        printf("\nNot enough memory for c[%d]\n", n);
        exit(1);
    }

    d = new(nothrow) double [n];
    if (!d) {
        printf("\n\nError in function polint():");
        printf("\nNot enough memory for d[%d]\n", n);
        exit(1);
    }

    for (i = 0; i < n; i++) {
        if ((dift = fabs(x - xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    *y = ya[ns--];
    for (m = 0; m < (n - 1); m++) {
        for (i = 0; i < n - m; i++) {
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            if ((den = ho - hp) < ZERO) {
                printf("\n\n Error in function polint(): ");
                printf("\nden = ho - hp = %4.1E -- too small\n", den);
                exit(1);
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        *y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
    delete [] d;
    delete [] c;

} // End: function polint()

/*
 ** The function
 **       rtbis()
 ** calculates a root between x1 and x2 of a function
 ** pointed to by (*func) using the method of bisection  
 ** The root is returned with an accuracy of +- xacc.
 */

#define MAXIT 40                   // max iterations

double rtbis(double (*func)(double), double x1, double x2, double xacc) {
    int j;
    double dx, f, fmid, xmid, rtb;

    f = (*func)(x1);
    fmid = (*func)(x2);
    if (f * fmid >= 0.0) {
        printf("\n\nError in function rtbis():");
        printf("\nroot in function must be within");
        printf("x1 = %f  and x2 = %f\n", x1, x2);
        exit(1);
    }
    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    for (j = 0; j < MAXIT; j++) {
        fmid = (*func)(xmid = rtb + (dx *= 0.5));
        if (fmid <= 0.0) rtb = xmid;
        if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    printf("\n\nError in function rtbis():"); // should never reach this point
    printf("\nToo many bisections!!!\n");
    exit(1);
}
#undef MAXIT
// End: function rtbis()

/*
 ** The function
 **      rtsec()
 ** calculates a root between x1 and x2 of a function
 ** pointed to by (*func) using the secant method.
 ** The root is returned with an accuracy of +- xacc.
 */

#define MAXIT 30      // max iterations 

double rtsec(double (*func)(double), double x1, double x2, double xacc) {
    int j;
    double fl, f, dx, swap, xl, rts;

    fl = (*func)(x1);
    f = (*func)(x2);
    if (fabs(fl) < fabs(f)) {
        rts = x1;
        xl = x2;
        swap = fl;
        fl = f;
        f = swap;
    } else {
        xl = x1;
        rts = x2;
    }
    for (j = 0; j < MAXIT; j++) {
        dx = (xl - rts) * f / (f - fl);
        xl = rts;
        fl = f;
        rts += dx;
        f = (*func)(rts);
        if (fabs(dx) < xacc || f == 0.0) return rts;
    }
    printf("\n\nError in function rtsec():"); // should never reach this point 
    printf("\nToo many iterations!!!\n");
    exit(1);
}
#undef MAXIT
// End: function rtsec()

/*
 ** The function
 **       rtnewt()
 ** calculates a root between x1 and x2 of a function pointed to
 ** by (*funcd) using the Newton-Raphson method. The user-defined
 ** function funcd() returns both the function value and its first
 ** derivative at the point x,
 ** The root is returned with an accuracy of +- xacc.
 */

#define MAXIT 20                       // max iterations

double rtnewt(void (*funcd)(double, double *, double *), double x1, double x2,
        double xacc) {
    int j;
    double df, dx, f, rtn;

    rtn = 0.5 * (x1 + x2); // initial guess 
    for (j = 0; j < MAXIT; j++) {
        (*funcd)(rtn, &f, &df);
        dx = f / df;
        rtn -= dx;
        if ((x1 - rtn) * (rtn - x2) < 0.0) {
            printf("\n\nError in function rtnewt():");
            printf("\nJump out of interval bracket\n");
            exit(1);
        }
        if (fabs(dx) < xacc) return rtn;
    }
    printf("\n\nError in function rtnewt():"); // should never reach this point
    printf("\nToo many iterations!!!\n");
    exit(1);
}
#undef MAXIT
// End: function rtnewt()

/*
 ** The function
 **       zbrent()
 ** calculates a root between x1 and x2 of a function
 ** pointed to by (*funcd) using the Brent's method.
 ** The root is returned with an accuracy of +- xacc.
 */

#define MAXIT 20                       // max iterations
#define EPS 3.0E-8

double zbrent(double (*func)(double), double x1, double x2, double xacc) {
    int iter;
    double a = x1, b = x2, c = x2, d, e, min1, min2;
    double fa = (*func)(a), fb = (*func)(b), fc, p, q, r, s, xacc1, xm;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        printf("\n\nError in function zbrent();");
        printf("\nRoot must be bracketed in zbrent()\n");
        exit(1);
    }

    fc = fb;
    for (iter = 1; iter <= MAXIT; iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (fabs(fc) < fabs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        xacc1 = 2.0 * EPS * fabs(b) + 0.5 * xacc;
        xm = 0.5 * (c - b);
        if (fabs(xm) <= xacc1 || fb == 0.0) return b;
        if (fabs(e) >= xacc1 && fabs(fa) > fabs(fb)) {
            s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) *(s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = fabs(p);
            min1 = 3.0 * xm * q - fabs(xacc1 * q);
            min2 = fabs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (fabs(d) > xacc1) b += d;
        else b += SIGN(xacc1, xm);
        fb = (*func)(b);
    }
    printf("\n\nError in function zbrent():"); // should never reach this point
    printf("\nToo many iterations!!!\n");
    exit(1);
}
#undef MAXIT
#undef EPS
// End: function zbrent()

/*
 ** The function
 **           ran0()
 ** is an "Minimal" random number generator of Park and Miller
 ** (see Numerical recipe page 279). Set or reset the input value
 ** idum to any integer value (except the unlikely value MASK)
 ** to initialize the sequence; idum must not be altered between
 ** calls for sucessive deviates in a sequence.
 ** The function returns a uniform deviate between 0.0 and 1.0.
 */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum) {
    long k;
    double ans;

    *idum ^= MASK;
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR*k;
    if (*idum < 0) *idum += IM;
    ans = AM * (*idum);
    *idum ^= MASK;
    return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

// End: function ran0() 

/*
 ** The function
 **           ran1()
 ** is an "Minimal" random number generator of Park and Miller
 ** (see Numerical recipe page 280) with Bays-Durham shuffle and
 ** added safeguards. Call with idum a negative integer to initialize;
 ** thereafter, do not alter idum between sucessive deviates in a
 ** sequence. RNMX should approximate the largest floating point value
 ** that is less than 1.
 ** The function returns a uniform deviate between 0.0 and 1.0
 ** (exclusive of end-point values).
 */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum) {
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR*k;
    if (*idum < 0) *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran1()

/*
 ** The function 
 **         ran2()
 ** is a long periode (> 2 x 10^18) random number generator of 
 ** L'Ecuyer and Bays-Durham shuffle and added safeguards.
 ** Call with idum a negative integer to initialize; thereafter,
 ** do not alter idum between sucessive deviates in a
 ** sequence. RNMX should approximate the largest floating point value
 ** that is less than 1.
 ** The function returns a uniform deviate between 0.0 and 1.0
 ** (exclusive of end-point values).
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum) {
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k*IR1;
    if (*idum < 0) *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran2()

/*
 ** The function
 **        ran3()
 ** returns a uniform random number deviate between 0.0 and 1.0. Set
 ** the idum to any negative value to initialize or reinitialize the
 ** sequence. Any large MBIG, and any small (but still large) MSEED
 ** can be substituted for the present values. 
 */

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum) {
    static int inext, inextp;
    static long ma[56]; // value 56 is special, do not modify
    static int iff = 0;
    long mj, mk;
    int i, ii, k;

    if (*idum < 0 || iff == 0) { // initialization
        iff = 1;

        mj = MSEED - (*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55] = mj; // initialize ma[55] 

        for (i = 1, mk = 1; i <= 54; i++) { // initialize rest of table 
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) mk += MBIG;
            mj = ma[ii];
        }

        for (k = 1; k <= 4; k++) { // randimize by "warming up" the generator
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        }

        inext = 0; // prepare indices for first generator number
        inextp = 31; // 31 is special
        *idum = 1;
    }

    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;
    mj = ma[inext] - ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext] = mj;
    return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

// End: function ran3()

double trapezoidal_rule(double a, double b, int n, double (*func)(double)) {
    double trapez_sum;
    double fa, fb, x, step;
    int j;
    step = (b - a) / ((double) n);
    fa = (*func)(a) / 2.;
    fb = (*func)(b) / 2.;
    trapez_sum = 0.;
    for (j = 1; j <= n - 1; j++) {
        x = j * step + a;
        trapez_sum += (*func)(x);
    }
    trapez_sum = (trapez_sum + fb + fa) * step;
    return trapez_sum;
} // end trapezoidal_rule 

double rectangle_rule(double a, double b, int n, double (*func)(double)) {
    double rectangle_sum;
    double fa, fb, x, step;
    int j;
    step = (b - a) / ((double) n);
    rectangle_sum = 0.;
    for (j = 0; j <= n; j++) {
        x = (j + 0.5) * step; // midpoint of a given rectangle
        rectangle_sum += (*func)(x); //  add value of function.
    }
    rectangle_sum *= step; //  multiply with step length.
    return rectangle_sum;
} // end rectangle_rule 

// Calculates the factorial of n

int factorial(int n) {
    if (n <= 1)
        return 1;
    else
        return n * factorial(n - 1);
}

/*
 ** The function
 **        ran3()
 ** returns a uniform random number deviate between 0.0 and 1.0. Set
 ** the idum to any negative value to initialize or reinitialize the
 ** sequence. Any large MBIG, and any small (but still large) MSEED
 ** can be substituted for the present values. 
 */
#if 1

void gauher(double x[], double w[], int n) {
    switch (n) {

        case 2:
            x[0] = -7.07106781186547524382e-01;
            x[1] = 7.07106781186547524382e-01;

            w[0] = 8.86226925452758013655e-01;
            w[1] = 8.86226925452758013655e-01;
            break;
        case 3:
            x[0] = -1.22474487139158904909e+00;
            x[1] = 0;
            x[2] = 1.22474487139158904909e+00;

            w[0] = 2.95408975150919337894e-01;
            w[1] = 1.18163590060367735158e+00;
            w[2] = 2.95408975150919337894e-01;
            break;
        case 4:
            x[0] = -1.65068012388578455585e+00;
            x[1] = -5.24647623275290317900e-01;
            x[2] = 5.24647623275290317900e-01;
            x[3] = 1.65068012388578455585e+00;

            w[0] = 8.13128354472451771398e-02;
            w[1] = 8.04914090005512836482e-01;
            w[2] = 8.04914090005512836482e-01;
            w[3] = 8.13128354472451771398e-02;
            break;
        case 20:
            x[0] = -5.38748089001123286199e+00;
            x[1] = -4.60368244955074427298e+00;
            x[2] = -3.94476404011562521040e+00;
            x[3] = -3.34785456738321632688e+00;
            x[4] = -2.78880605842813048055e+00;
            x[5] = -2.25497400208927552311e+00;
            x[6] = -1.73853771211658620678e+00;
            x[7] = -1.23407621539532300786e+00;
            x[8] = -7.37473728545394358719e-01;
            x[9] = -2.45340708300901249903e-01;
            x[10] = 2.45340708300901249903e-01;
            x[11] = 7.37473728545394358719e-01;
            x[12] = 1.23407621539532300786e+00;
            x[13] = 1.73853771211658620678e+00;
            x[14] = 2.25497400208927552311e+00;
            x[15] = 2.78880605842813048055e+00;
            x[16] = 3.34785456738321632688e+00;
            x[17] = 3.94476404011562521040e+00;
            x[18] = 4.60368244955074427298e+00;
            x[19] = 5.38748089001123286199e+00;

            w[0] = 2.22939364553415129254e-13;
            w[1] = 4.39934099227318055366e-10;
            w[2] = 1.08606937076928169398e-07;
            w[3] = 7.80255647853206369398e-06;
            w[4] = 2.28338636016353967260e-04;
            w[5] = 3.24377334223786183217e-03;
            w[6] = 2.48105208874636108814e-02;
            w[7] = 1.09017206020023320014e-01;
            w[8] = 2.86675505362834129720e-01;
            w[9] = 4.62243669600610089640e-01;
            w[10] = 4.62243669600610089640e-01;
            w[11] = 2.86675505362834129720e-01;
            w[12] = 1.09017206020023320014e-01;
            w[13] = 2.48105208874636108814e-02;
            w[14] = 3.24377334223786183217e-03;
            w[15] = 2.28338636016353967260e-04;
            w[16] = 7.80255647853206369398e-06;
            w[17] = 1.08606937076928169398e-07;
            w[18] = 4.39934099227318055366e-10;
            w[19] = 2.22939364553415129254e-13;
            break;
        case 127:
            x[0] = -15.22833814816731;
            x[1] = -14.66959515883396;
            x[2] = -14.20908599528489;
            x[3] = -13.79972229021166;
            x[4] = -13.42351859007094;
            x[5] = -13.07120866047461;
            x[6] = -12.73723565241570;
            x[7] = -12.41793937886972;
            x[8] = -12.11074902094773;
            x[9] = -11.81377219826773;
            x[10] = -11.52556511257270;
            x[11] = -11.24499458378554;
            x[12] = -10.97115056984025;
            x[13] = -10.70328820102748;
            x[14] = -10.44078795777277;
            x[15] = -10.18312747345035;
            x[16] = -9.929861049511418;
            x[17] = -9.680604441247482;
            x[18] = -9.435023338988188;
            x[19] = -9.192824498846019;
            x[20] = -8.953748810856546;
            x[21] = -8.717565808707640;
            x[22] = -8.484069268983252;
            x[23] = -8.253073645445717;
            x[24] = -8.024411151470339;
            x[25] = -7.797929351387017;
            x[26] = -7.573489155608352;
            x[27] = -7.350963139226906;
            x[28] = -7.130234122035071;
            x[29] = -6.911193961546574;
            x[30] = -6.693742520875827;
            x[31] = -6.477786781164532;
            x[32] = -6.263240074273735;
            x[33] = -6.050021416141988;
            x[34] = -5.838054924877421;
            x[35] = -5.627269310546483;
            x[36] = -5.417597425924316;
            x[37] = -5.208975869315402;
            x[38] = -5.001344632038633;
            x[39] = -4.794646784376496;
            x[40] = -4.588828194769841;
            x[41] = -4.383837277846472;
            x[42] = -4.179624767535203;
            x[43] = -3.976143512067340;
            x[44] = -3.773348288125052;
            x[45] = -3.571195631778216;
            x[46] = -3.369643684171737;
            x[47] = -3.168652050195365;
            x[48] = -2.968181668595591;
            x[49] = -2.768194692182406;
            x[50] = -2.568654376947349;
            x[51] = -2.369524979049037;
            x[52] = -2.170771658741149;
            x[53] = -1.972360390419503;
            x[54] = -1.774257878051679;
            x[55] = -1.576431475326778;
            x[56] = -1.378849109926179;
            x[57] = -1.181479211370068;
            x[58] = -0.9842906419402724;
            x[59] = -0.7872526302182505;
            x[60] = -0.5903347068094211;
            x[61] = -0.3935066418513014;
            x[62] = -0.1967383839242327;
            x[63] = -0.1140913802913524E-16;
            x[64] = 0.1967383839242332;
            x[65] = 0.3935066418513009;
            x[66] = 0.5903347068094202;
            x[67] = 0.7872526302182506;
            x[68] = 0.9842906419402723;
            x[69] = 1.181479211370069;
            x[70] = 1.378849109926177;
            x[71] = 1.576431475326779;
            x[72] = 1.774257878051679;
            x[73] = 1.972360390419502;
            x[74] = 2.170771658741150;
            x[75] = 2.369524979049042;
            x[76] = 2.568654376947349;
            x[77] = 2.768194692182402;
            x[78] = 2.968181668595590;
            x[79] = 3.168652050195361;
            x[80] = 3.369643684171739;
            x[81] = 3.571195631778221;
            x[82] = 3.773348288125053;
            x[83] = 3.976143512067337;
            x[84] = 4.179624767535202;
            x[85] = 4.383837277846474;
            x[86] = 4.588828194769836;
            x[87] = 4.794646784376495;
            x[88] = 5.001344632038628;
            x[89] = 5.208975869315393;
            x[90] = 5.417597425924324;
            x[91] = 5.627269310546481;
            x[92] = 5.838054924877422;
            x[93] = 6.050021416141982;
            x[94] = 6.263240074273729;
            x[95] = 6.477786781164535;
            x[96] = 6.693742520875829;
            x[97] = 6.911193961546577;
            x[98] = 7.130234122035069;
            x[99] = 7.350963139226908;
            x[100] = 7.573489155608353;
            x[101] = 7.797929351387015;
            x[102] = 8.024411151470339;
            x[103] = 8.253073645445719;
            x[104] = 8.484069268983250;
            x[105] = 8.717565808707629;
            x[106] = 8.953748810856544;
            x[107] = 9.192824498846035;
            x[108] = 9.435023338988174;
            x[109] = 9.680604441247466;
            x[110] = 9.929861049511425;
            x[111] = 10.18312747345035;
            x[112] = 10.44078795777277;
            x[113] = 10.70328820102749;
            x[114] = 10.97115056984024;
            x[115] = 11.24499458378554;
            x[116] = 11.52556511257271;
            x[117] = 11.81377219826772;
            x[118] = 12.11074902094773;
            x[119] = 12.41793937886973;
            x[120] = 12.73723565241570;
            x[121] = 13.07120866047460;
            x[122] = 13.42351859007094;
            x[123] = 13.79972229021168;
            x[124] = 14.20908599528486;
            x[125] = 14.66959515883398;
            x[126] = 15.22833814816733;
            w[0] = 0.1250449757709384E-100;
            w[1] = 0.1727279805947516E-93;
            w[2] = 0.8932168157226119E-88;
            w[3] = 0.7730618524115457E-83;
            w[4] = 0.2014395765271519E-78;
            w[5] = 0.2150371473367845E-74;
            w[6] = 0.1134192420862984E-70;
            w[7] = 0.3348913901189688E-67;
            w[8] = 0.6048654896421865E-64;
            w[9] = 0.7137509294658209E-61;
            w[10] = 0.5788456337506695E-58;
            w[11] = 0.3358116622396395E-55;
            w[12] = 0.1439464194929904E-52;
            w[13] = 0.4682180838336620E-50;
            w[14] = 0.1181705444072102E-47;
            w[15] = 0.2358165915608093E-45;
            w[16] = 0.3781442794092380E-43;
            w[17] = 0.4941103111592451E-41;
            w[18] = 0.5325530377558451E-39;
            w[19] = 0.4785439068028819E-37;
            w[20] = 0.3619188344606345E-35;
            w[21] = 0.2323208338641435E-33;
            w[22] = 0.1275333141104785E-31;
            w[23] = 0.6027775385095142E-30;
            w[24] = 0.2467977324185458E-28;
            w[25] = 0.8801956769196604E-27;
            w[26] = 0.2748248921212509E-25;
            w[27] = 0.7546821890331683E-24;
            w[28] = 0.1830313463633893E-22;
            w[29] = 0.3935599086098561E-21;
            w[30] = 0.7529316163881436E-20;
            w[31] = 0.1285799778676248E-18;
            w[32] = 0.1965932688850777E-17;
            w[33] = 0.2698651190729676E-16;
            w[34] = 0.3334441430330234E-15;
            w[35] = 0.3717330312526506E-14;
            w[36] = 0.3747395447295843E-13;
            w[37] = 0.3423009449350244E-12;
            w[38] = 0.2838530372508185E-11;
            w[39] = 0.2140692029052080E-10;
            w[40] = 0.1470633127347631E-09;
            w[41] = 0.9217394096771948E-09;
            w[42] = 0.5278166393713740E-08;
            w[43] = 0.2765049704503656E-07;
            w[44] = 0.1326785584258085E-06;
            w[45] = 0.5838094427629430E-06;
            w[46] = 0.2358156172484939E-05;
            w[47] = 0.8752446803455428E-05;
            w[48] = 0.2987679053600162E-04;
            w[49] = 0.9387443572036433E-04;
            w[50] = 0.2717076262801515E-03;
            w[51] = 0.7249392974272461E-03;
            w[52] = 0.1784120832681892E-02;
            w[53] = 0.4052485518617183E-02;
            w[54] = 0.8500026304180853E-02;
            w[55] = 0.1647114224166119E-01;
            w[56] = 0.2949929624830547E-01;
            w[57] = 0.4884738711445193E-01;
            w[58] = 0.7480798976881646E-01;
            w[59] = 0.1059852050812395;
            w[60] = 0.1389394530909479;
            w[61] = 0.1685623607426031;
            w[62] = 0.1892784958017932;
            w[63] = 0.1967334068888451;
            w[64] = 0.1892784958017945;
            w[65] = 0.1685623607426027;
            w[66] = 0.1389394530909484;
            w[67] = 0.1059852050812391;
            w[68] = 0.7480798976881628E-01;
            w[69] = 0.4884738711445169E-01;
            w[70] = 0.2949929624830535E-01;
            w[71] = 0.1647114224166098E-01;
            w[72] = 0.8500026304180756E-02;
            w[73] = 0.4052485518617234E-02;
            w[74] = 0.1784120832681892E-02;
            w[75] = 0.7249392974272203E-03;
            w[76] = 0.2717076262801586E-03;
            w[77] = 0.9387443572036826E-04;
            w[78] = 0.2987679053600191E-04;
            w[79] = 0.8752446803455300E-05;
            w[80] = 0.2358156172484923E-05;
            w[81] = 0.5838094427629397E-06;
            w[82] = 0.1326785584258084E-06;
            w[83] = 0.2765049704503730E-07;
            w[84] = 0.5278166393713601E-08;
            w[85] = 0.9217394096772193E-09;
            w[86] = 0.1470633127347794E-09;
            w[87] = 0.2140692029052102E-10;
            w[88] = 0.2838530372508322E-11;
            w[89] = 0.3423009449350289E-12;
            w[90] = 0.3747395447295533E-13;
            w[91] = 0.3717330312526600E-14;
            w[92] = 0.3334441430330258E-15;
            w[93] = 0.2698651190729916E-16;
            w[94] = 0.1965932688850665E-17;
            w[95] = 0.1285799778676250E-18;
            w[96] = 0.7529316163881871E-20;
            w[97] = 0.3935599086098257E-21;
            w[98] = 0.1830313463633751E-22;
            w[99] = 0.7546821890331568E-24;
            w[100] = 0.2748248921212469E-25;
            w[101] = 0.8801956769197440E-27;
            w[102] = 0.2467977324185710E-28;
            w[103] = 0.6027775385094256E-30;
            w[104] = 0.1275333141104785E-31;
            w[105] = 0.2323208338641707E-33;
            w[106] = 0.3619188344606339E-35;
            w[107] = 0.4785439068028225E-37;
            w[108] = 0.5325530377558222E-39;
            w[109] = 0.4941103111593767E-41;
            w[110] = 0.3781442794091814E-43;
            w[111] = 0.2358165915608099E-45;
            w[112] = 0.1181705444072005E-47;
            w[113] = 0.4682180838336191E-50;
            w[114] = 0.1439464194929855E-52;
            w[115] = 0.3358116622396488E-55;
            w[116] = 0.5788456337505780E-58;
            w[117] = 0.7137509294657940E-61;
            w[118] = 0.6048654896421083E-64;
            w[119] = 0.3348913901189406E-67;
            w[120] = 0.1134192420862977E-70;
            w[121] = 0.2150371473368259E-74;
            w[122] = 0.2014395765271411E-78;
            w[123] = 0.7730618524113554E-83;
            w[124] = 0.8932168157229719E-88;
            w[125] = 0.1727279805947236E-93;
            w[126] = 0.1250449757709172E-100;
            break;
    }
}
#endif

bool fexists(const char *filename) {
    ifstream ifile(filename);
    return ifile;
}

#if 0
void gauher(double *x, double *w, int n) {
    const double eps = 1.0e-17;
    const double pim4 = 0.7511255444649425;
    const int maxit = 10;
    int i, its, j, m;
    double p1, p2, p3, pp, z, z1;

    m = (n + 1) / 2;
    for (i = 0; i < m; i++) {
        if (i == 0) {
            z = sqrt(double(2 * n + 1)) - 1.85575 * pow(double(2 * n + 1), -0.16667);
        } else if (i == 1) {
            z -= 1.14 * pow(double(n), 0.426) / z;
        } else if (i == 2) {
            z = 1.86 * z - 0.86 * x[0];
        } else if (i == 3) {
            z = 1.91 * z - 0.91 * x[1];
        } else {
            z = 2.0 * z - x[i - 2];
        }
        for (its = 0; its < maxit; its++) {
            p1 = pim4;
            p2 = 0.0;
            for (j = 0; j < n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0 / (j + 1)) * p2 - sqrt(double(j) / (j + 1)) * p3;
            }
            pp = sqrt(double(2 * n)) * p2;
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= eps) break;
        }
        if (its >= maxit) {
            cout << "too many iterations in gauher";
            exit(1);
        }
        x[i] = z;
        x[n - 1 - i] = -z;
        w[i] = 2.0 / (pp * pp);
        w[n - 1 - i] = w[i];
    }
  
    /*
    //    void nrerror(char error_text[]);
    int maxit = 10;
    double eps = 3.0e-14;
    int i, its, j, m;
    double p1, p2, p3, pp, z, z1;

    m = (n + 1) / 2;
    for (i = 1; i <= m; i++) {
        if (i == 1) {
            z = sqrt((double) (2 * n + 1)) - 1.85575 * pow((double) (2 * n + 1), -0.16667);
        } else if (i == 2) {
            z -= 1.14 * pow((double) n, 0.426) / z;
        } else if (i == 3) {
            z = 1.86 * z - 0.86 * x[0];
            //            z = 1.86 * z - 0.86 * x[1];
        } else if (i == 4) {
            z = 1.91 * z - 0.91 * x[1];
            //            z = 1.91 * z - 0.91 * x[2];
        } else {
            z = 2.0 * z - x[i - 1];
            //            z = 2.0 * z - x[i - 2];
        }
        for (its = 1; its <= maxit; its++) {
            p1 = PIM4;
            p2 = 0.0;
            for (j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0 / j) * p2 - sqrt(((double) (j - 1)) / j) * p3;
            }
            pp = sqrt((double) 2 * n) * p2;
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= eps) break;
        }
        if (its > maxit) {
            cout << "too many iterations in gauher";
            exit(1);
        }
        //        if (its > MAXIT) nrerror("too many iterations in gauher");
        x[i - 1] = z;
        x[n - i] = -z;
        w[i - 1] = 2.0 / (pp * pp);
        w[n - i] = w[i - 1];

        //        x[i] = z;
        //        x[n + 1 - i] = -z;
        //        w[i] = 2.0 / (pp * pp);
        //        w[n + 1 - i] = w[i];
    }
     */
}
#endif