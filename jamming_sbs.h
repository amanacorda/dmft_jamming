//
//  jamming.h
//  
//
//  Created by Alessandro Manacorda on 5/20/20.
//

#ifndef active_h
#define active_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "random.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#define FLEN 50
#define INI_CHECK 630

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
    double y ;
    double h ;
    double f ;
    double eta ;
    double mu ;
    double v[3] ;
} traj_t ;

typedef struct {
    double G ;
    double k ;
    double chiR ;
    double chiMR ;
    double E ;
    double P ;
    double err ;
} vec_t ;

typedef struct {
    double R ;
    double C ;
    double MR ;
    double MC ;
    double Minv ;
    double Vn ;
    double MCtest ;
    double Gtest ;
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
    double G0 ;
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
    int N ;
} phys_t ;

typedef struct {
    double fc1 ;
    double fc2 ;
    double fc3 ;
    double dtaup ;
    double dt ;
    double dt2 ;
    double sM0 ;
    double eps ;
    double Teff ;
    int N ;
    int N1 ;
    int m0 ;
} coeffs_t ;

typedef struct {
    double k1 ;
    double *MR1 ;
    double *MC1 ;
    double E1 ;
    double P1 ;
} increments_t ;

typedef struct {
    traj_t **X ;
    double **H ;
    vec_t *v ;
    mat_t *M ;
    hweights_t *hw ;
} dyn_t ;

void init_par(int,char**,phys_t*,hint_t*,ini_t*) ;
void init_coeffs(phys_t,ini_t,coeffs_t*) ;
void init_dyn(dyn_t*,phys_t,hint_t*,ini_t) ;
void init_traj(dyn_t*,phys_t,hint_t,double,int) ;
void CR_step(vec_t*,mat_t*,coeffs_t,int) ;
void h0int(dyn_t*,hint_t,coeffs_t*,int,int,int) ;
void traj_step(dyn_t*,increments_t*,hint_t,coeffs_t,double,int,int,int) ;
//void check(traj);
void resize(dyn_t*,hint_t*,int) ;
double det(double**,int) ;
void matrix_inversion(mat_t*,int) ;
void block_matrix_inversion(mat_t*,int,int) ;
double check_id(mat_t*,int,int) ;
void stabilize(mat_t*,int) ;
void Cholesky_inversion(mat_t*,int) ;
void LU_inversion(mat_t*,int) ;
double check_id_gsl(gsl_matrix*,gsl_matrix*,int) ;
void stabilize_cholesky(mat_t*,int) ;
//void print_matrix(double*,char*,double,int) ;
void norm_test(dyn_t*,hint_t,int) ;
void print_output(vec_t*,mat_t*,phys_t) ;
void filename(char*,char*,phys_t) ;
void print_info(phys_t,hint_t,ini_t) ;
void free_dyn(dyn_t*,int,int) ;
void print_coeffs(coeffs_t) ;


#endif /* active_h */
