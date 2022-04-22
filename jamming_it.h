//
//  jamming_it.h
//  
//
//  Created by Alessandro Manacorda on 1/27/21.
//

#ifndef jamming_it_h
#define jamming_it_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> 
#include "random.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <omp.h>

#define FLEN 50
#define INI_CHECK 630

#ifndef N
#define N 301
#endif
#define N1 (int) (N*(N+1)/2)

#define CHECKS 1
#define CHECK_MATRICES 1
#define STABILIZATION 1

static inline double V(double h, int n) {
    return (h*h/2.*(n==0)+h*(n==1)+(n==2))*(h<0) ;
}

static inline double bV0(double h, double b) {
    if(h>=0) return 0. ;
    else if(isinf(b)) return INFINITY ;
    else return (b*V(h,0)) ;
}

static inline int tri(int n, int m) {
    if(n>=m) return (int)(n*(n+1)/2+m) ;
    else return tri(m,n) ;
}

static inline double rel_err(double Fn, double Ft) {
    if(Ft!=0) return fabs(Fn/Ft-1.) ;
    else return fabs(Ft/Fn-1.) ;
}


typedef struct {
    double y[N] ;
    double h[N] ;
    //double f ;
    double eta[N] ;
    //double mu ;
    double v[N][3] ;
    double F[N1] ;
} traj_t ;

typedef struct {
    //double k[N] ;
    //double chiR[N] ;
    double E[N] ;
    double P[N] ;
    double err[N] ;
} vec_t ;

typedef struct {
    double R[N1] ;
    double C[N1] ;
    //double MR[N1] ;
    double chi[N1] ;
    double MC[N1] ;
    //double Minv ;
    //double Vn ;
    double MCtest[N1] ;
} mat_t ;

typedef struct {
    double hmin ;
    double hmax ;
    int NH ;
    double dh0 ;
    double z ;
    double norm ;
    int NP ;
    double threshold ;
    double normNR ;
} hint_t ;

typedef struct {
    double h0 ;
    double wh0 ;
} hweights_t ;

typedef struct {
    double k0 ;
    double MR0 ;
    double MC0 ;
    double E0 ;
    double P0 ;
} ini_t ;

typedef struct {
    double phi ;
    double eps ;
    double Teff ;
    double taup ;
    double b0 ;
    double dt ;
    double alpha ;
    //int N ;
} phys_t ;

typedef struct {
    double dt ;
    double dt2 ;
    double sM0 ;
    double eps ;
    double Teff ;
    int m0 ;
    int noise ;
    double noise_threshold ;
    double noise_mult ;

} coeffs_t ;

typedef struct {
    //double k1[N] ;
    //double MR1[N1] ;
    double chi1[N1] ;
    double MC1[N1] ;
    double E1[N] ;
    double P1[N] ;
} increments_t ;

typedef struct {
    traj_t X ;
    //double H ;
    vec_t v ;
    mat_t M ;
    hweights_t *hw ;
} dyn_t ;

void init_dyn_it(dyn_t*,phys_t,hint_t*,ini_t) ;
void init_par(int,char**,phys_t*,hint_t*,ini_t*) ;
void init_coeffs(phys_t,ini_t,coeffs_t*) ;
void print_info(phys_t,hint_t,ini_t) ;

void eigensolve(mat_t,gsl_vector*,gsl_matrix*);
void noise_generation(gsl_vector*,gsl_vector*,gsl_matrix*);
void noise_test(dyn_t*) ;
void noise_norm(mat_t*,double);

//void init_traj_it(dyn_t*,phys_t,hint_t,int) ;
void init_temp(increments_t*) ;
void measure_kernels(dyn_t*,increments_t*,double) ;
void update_kernels(dyn_t*,increments_t*,double*,phys_t) ;

void CR_solve(vec_t*,mat_t*,coeffs_t) ;
void traj_solve(dyn_t*,coeffs_t,double) ;

void print_output(vec_t,mat_t,phys_t,char*) ;
void filename(char*,char*,phys_t,char*) ;

#endif /* jamming_it_h */
