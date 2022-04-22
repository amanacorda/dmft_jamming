//
//  jamming_sbs.c
//  
//
//  Created by Alessandro Manacorda on 5/19/20.
//


#include "jamming_sbs.h"

int main(int argc, char **argv) {
    
    phys_t p0 ; hint_t hi0 ; ini_t i0 ;         // Declaration of physical parameters, integration constants, initial conditions
    init_par(argc-1,argv,&p0,&hi0,&i0) ;        // Initialization of parameters
    coeffs_t c0 ;                               // Declaration of numerical coefficients
    init_coeffs(p0,i0,&c0) ;                    // Initialization
    
    dyn_t *D = malloc(sizeof(dyn_t)) ;          // Declaration of dynamical variables
    init_dyn(D,p0,&hi0,i0) ;                    // Initialization
    int n, N=p0.N, N1=c0.N1, freq=10 ;
    
    FILE *fh ; fh = fopen("traj.dat","w") ;
    fprintf(fh,"%lf",0.);
    for(int nh=0;nh<hi0.NH;nh++) fprintf(fh,"\t%21.15e",D->X[0][nh*hi0.NP].h) ;
    fprintf(fh,"\n");
    
    clock_t start = clock() ;                   // Initial time
    for(n=0;n<N-1;n++) {
        CR_step(D->v,D->M,c0,n);                // Evolution of C and R
        
        h0int(D,hi0,&c0,n+1,n,0) ;                // Integration over h0 and realizations
        fprintf(fh,"%lf",(n+1)*p0.dt);
        for(int nh=0;nh<hi0.NH;nh++) fprintf(fh,"\t%21.15e",D->X[n+1][nh*hi0.NP].h) ;
        fprintf(fh,"\n");
        if((n+1)%freq==0) {                     // Prints
            clock_t end = clock() ;
            printf("\nTime to execute %d steps = %f seconds\n\n",n+1,(float)(end-start)/CLOCKS_PER_SEC) ;
            fflush(stdout);
            print_output(D->v,D->M,p0) ;        // Output print
        }
       
        block_matrix_inversion(D->M,n,c0.m0) ;
    }
    
    printf("Fatto\n");
    
    //norm_test(D,hi0,N) ;                // Normalization of correlations
    
    print_output(D->v,D->M,p0) ;        // Output print
    fclose(fh) ;
    free_dyn(D,N,N1) ;                  // Free the dyn memory
    
    return 0 ;
    
}



// FUNCTION DEFINITIONS //


// Initialization of physical parameters, integration constants and initial conditions : from line argument or from code

void init_par(int a, char **argv, phys_t *par, hint_t *hi, ini_t *I0) {
    if(a)
    {
        if(a==9)
        {
            par->phi = atof(argv[1]) ;          // Rescaled packing fraction
            par->eps = atof(argv[2]) ;          // Potential energy scale
            //par->Teff = atof(argv[3]) ;         // Effective temperature
            //par->taup = atof(argv[4]) ;         // Persistence time
            par->b0 = atof(argv[3]) ;           // Initial beta
            par->dt = atof(argv[4]) ;           // Time step
            par->N = atoi(argv[5]) ;            // Number of time steps
            hi->dh0 = atof(argv[6]) ;           // Integration step
            hi->z = atof(argv[7]) ;             // Rescaled limit of integration
            hi->threshold = atof(argv[8]) ;    // Error threshold
            hi->NP = atoi(argv[9]) ;           // Number of realizations
        }
        else
        {
            printf("Wrong number of arguments!!\n") ;
            exit(EXIT_FAILURE) ;
        }
    }
    else {
        par->phi = 1.0 ;
        par->eps = 10. ;
        //par->Teff = 1.   ;
        //par->taup = 1. ;
        par->b0 = 0. ;
        par->dt = 0.05 ;
        par->N = 301 ;
        hi->dh0 = 0.01 ;
        hi->z = 50. ;
        hi->threshold = 1e-3 ;
        hi->NP = 100 ;
    }
    //double Tr=1./(par->eps*par->b0) ;
    //double wr=exp(Tr/2.)*erfc(sqrt(Tr/2.)) ;
    //printf("Tr = %lf\t wr = %lf\n",Tr,wr) ;
    /*
    I0->G0 = par->Teff/par->taup ;
    if(isfinite(par->b0)) I0->k0 = (par->phi/2.)*(sqrt(M_PI/(2.*Tr))*(1.+Tr)*wr - 1.)/par->b0 ;
    else if isinf(par->b0) I0->k0 = 0. ;
    I0->MR0 = (par->phi/2.)*par->eps*par->eps*sqrt(M_PI*Tr/2.)*wr ;
    I0->MC0 = I0->k0/par->b0 ;
    I0->E0 = Tr*I0->k0/2. ;
    I0->P0 = (I0->MR0/par->eps - I0->k0)*par->Teff ; // CHECK! What happens with beta??
     */
    I0->G0 = 0. ;
    if(par->b0==0){
        I0->k0 = 0. ;
        I0->MR0 = par->phi*par->eps*par->eps/2. ;
        I0->MC0 = 2.*I0->MR0 ;
        I0->E0 = par->phi*par->eps/2. ;
        I0->P0 = I0->E0 ;
    }
    else if(par->b0>0){
        double x = 1./(2.*par->b0*par->eps) ;
        double y = exp(x)*sqrt(M_PI*x)*erfc(sqrt(x)) ;
        I0->k0 = par->phi*(-1.+(1.+1./(2.*x))*y)/(2.*par->b0) ;
        I0->MR0 = par->phi*par->eps*par->eps*y/2. ;
        I0->MC0 = I0->k0/par->b0 ;
        I0->E0 = x*I0->k0 ;
        I0->P0 = par->phi*(1-y)/(2.*par->b0) ;
    }
    else if(par->b0<0){
        printf("Warning: negative temperature?!?\n\n");
        exit(EXIT_FAILURE) ;
    }
    
    //hi->hmin = -5.*sqrt(Tr) ;
    //hi->hmax = 2.*hi->z*par->dt/sqrt(I0->G0+I0->MC0) ;
    hi->hmin = -10. ;
    hi->hmax = sqrt(I0->MC0)*par->dt*hi->z + I0->MC0*par->dt*par->dt/2. ;
    //hi->hmax = sqrt(I0->G0+I0->MC0)*par->dt*(hi->z + sqrt(I0->G0+I0->MC0)/2.) ;
    hi->NH = (int)((hi->hmax-hi->hmin)/hi->dh0) ;
    hi->norm = par->phi*hi->dh0/(2.*hi->NP) ;
    hi->normNR = 1./((double)(hi->NH*hi->NP)) ;
}


// Initialization of numerical coefficients

void init_coeffs(phys_t par, ini_t I0, coeffs_t *c) {
    //c->dtaup = par.dt/par.taup ;
    c->dt = par.dt ;
    c->dt2 = par.dt*par.dt ;
    //c->fc1 = 1.-0.5*c->dtaup ;
    //c->fc2 = sqrt(2*par.Teff*par.dt)/par.taup ;
    //c->fc3 = 1./(1.+0.5*c->dtaup) ;
    c->sM0 = sqrt(I0.MC0) ;
    c->eps = par.eps ;
    //c->Teff = par.Teff ;
    c->N = par.N ;
    c->N1 = (int)(par.N*(par.N+1)/2) ;
    //c->m0 = isinf(par.b0) ;
    c->m0=0. ;
}


// Initialization of dynamical variables

void init_dyn(dyn_t *D, phys_t par, hint_t *hi, ini_t I0) {
    int NR = hi->NH*hi->NP, N=par.N ;
    int n, N1 = (int)(N*(N+1)/2) ;
    D->X = (traj_t**) calloc(N,sizeof(traj_t*)) ;
    for(n=0;n<N;n++) D->X[n] = (traj_t*) calloc(NR,sizeof(traj_t)) ;
    D->H = (double**) calloc(N1,sizeof(double*)) ;
    for(n=0;n<N1;n++) D->H[n] = (double*) calloc(NR,sizeof(double)) ;
    D->v = (vec_t*) calloc(N,sizeof(vec_t)) ;
    D->M = (mat_t*) calloc(N1,sizeof(mat_t)) ;
    D->hw = (hweights_t*) calloc(hi->NH,sizeof(hweights_t)) ;
    
    //double dtau = par.dt/par.taup ;
    for(int n=0;n<N;n++)
    {
        //D->v[n].G = I0.G0*exp(-n*dtau) ;
        D->v[n].k = 0. ;
        D->v[n].E = 0. ;
        D->v[n].P = 0. ;
    }
    for(int n=0;n<N1;n++) {
        D->M[n].R=0. ;
        D->M[n].C=0. ;
        D->M[n].MR=0. ;
        D->M[n].MC=0. ;
        D->M[n].Minv=0. ;
        D->M[n].Vn=0. ;
        D->M[n].MCtest = 0. ;
        D->M[n].Gtest = 0. ;
    }
    double sG0 = sqrt(I0.G0) ;
    D->M[0].R = 0.5 ;
    D->M[0].C = 0.;
    D->v[0].k = I0.k0 ;
    D->M[0].MR = I0.MR0 ;
    D->M[0].MC = I0.MC0 ;
    D->v[0].chiR = 0. ;
    D->v[0].chiMR = 0. ;
    D->M[0].Minv = 1./D->M[0].MC ;
    D->M[0].Vn = D->M[0].Minv ;
    D->v[0].E = I0.E0 ;
    D->v[0].P = I0.P0 ;
    
    //double Tr = 1./(par.eps*par.b0) ;
    double e0=0., p0=0., k0=0., MR0=0., MC0=0., test=0. ;
    // TEST OF THE INTEGRATION OVER h0
    FILE *fh ; fh = fopen("hw0.dat","w");
    for(int nh=0;nh<hi->NH;nh++) {
        init_traj(D,par,*hi,sG0,nh) ;
        double x ;
        fprintf(fh,"%lf\t%21.15e\t",D->hw[nh].h0,D->hw[nh].wh0);
        x = hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[0] ; fprintf(fh,"%21.15e\t",x);
        e0 += x ;
        x = -hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[1] ; fprintf(fh,"%21.15e\t",x);
        p0 += x ;
        x = hi->NP*D->hw[nh].wh0*(D->X[0][hi->NP*nh].v[2]+D->X[0][hi->NP*nh].v[1]) ; fprintf(fh,"%21.15e\t",x);
        k0 += x ;
        x = hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[2]*D->H[0][hi->NP*nh] ; fprintf(fh,"%21.15e\t",x);
        MR0 += x ;
        x = hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[1]*D->X[0][hi->NP*nh].v[1] ; fprintf(fh,"%21.15e\t",x);
        MC0 += x ;
        x = hi->NP*D->hw[nh].wh0 ;
        test += x ; fprintf(fh,"%21.15e\t",x);
        fprintf(fh,"\n");
        /*
        e0 += hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[0] ;
        p0 -= hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[1] ;
        k0 += hi->NP*D->hw[nh].wh0*(D->X[0][hi->NP*nh].v[2]+D->X[0][hi->NP*nh].v[1]) ;
        MR0 += hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[2]*D->H[0][hi->NP*nh] ;
        MC0 += hi->NP*D->hw[nh].wh0*D->X[0][hi->NP*nh].v[1]*D->X[0][hi->NP*nh].v[1] ;
         */
    }
    fclose(fh);
    printf("k0 = %lf\t MR0 = %lf\t MC0 = %lf\n",k0,MR0,MC0) ;
    printf("e0 = %lf\t p0 = %lf\n",e0,p0) ;
    printf("test = %lf\n",test);
    printf("k0 = %lf\t MR0 = %lf\t MC0 = %lf\n",I0.k0,I0.MR0,I0.MC0) ;
    printf("e0 = %lf\t p0 = %lf\n",I0.E0,I0.P0) ;
    printf("k_err = %lf\t MR_err = %lf\t MC_err = %lf\n",rel_err(k0,I0.k0),rel_err(MR0,I0.MR0),rel_err(MC0,I0.MC0)) ;
    printf("e_err = %lf\t p_err = %lf\n",rel_err(e0,I0.E0),rel_err(p0,I0.P0)) ;
    double c1=rel_err(k0,I0.k0) ;
    c1 = fmax(c1,rel_err(MR0,I0.MR0)) ;
    c1 = fmax(c1,rel_err(MC0,I0.MC0)) ;
    if(isnormal(c1)) hi->threshold = fmax(c1,hi->threshold) ;
    printf("th = %lf\n",hi->threshold) ;
    D->v[0].err = c1 ;
    print_info(par,*hi,I0) ;                     // Print file with informations
}

void init_traj(dyn_t *D, phys_t par, hint_t hi, double sG0, int nh) {
    double h0 = hi.hmin + (nh+0.5)*hi.dh0 ;
    D->hw[nh].h0 = h0 ;
    double vv[3] ;
    for(int m=0;m<3;m++) vv[m] = par.eps*V(h0,m) ;
    for(int np=0;np<hi.NP;np++) {
        int nr = hi.NP*nh+np ;
        D->X[0][nr].y = 0. ;
        D->X[0][nr].h = h0 ;
        //D->X[0][nr].f = 0. ;
        //D->X[0][nr].f = sG0*GAUSSGEN ; D->M[0].Gtest += (D->X[0][nr].f)*(D->X[0][nr].f) ;
        for(int m=0;m<3;m++) D->X[0][nr].v[m] = vv[m] ;
        D->H[0][nr] = D->X[0][nr].v[2] ;
    }
    D->hw[nh].wh0 = hi.norm*exp(h0-par.eps*bV0(h0,par.b0)) ;
    
    //D->hw[nh].wh0 = hi.norm*exp(h0) ;
    //printf("V(%lf) = %lf\t",h0,V(h0,0)) ;
    //printf("w[%lf] = %lf\n",h0,D->hw[nh].wh0) ;
}


// Euler algorithm for the evolution of C and R

void CR_step(vec_t *v, mat_t *M, coeffs_t c, int n) {
    
    int n0=tri(n,0), n1=tri(n+1,0), nn=tri(n+1,n+1), l, m ;
    for(m=0;m<=n;m++) {
        double Rint=0., C1int=0., C2int=0. ;
        for(l=m+1;l<=n;l++) {
            int lm = tri(l,m) ;
            Rint += M[n0+l].MR*M[lm].R ;
            C1int += M[n0+l].MR*M[lm].C ;
        }
        M[n1+m].R = M[n0+m].R*(1.-v[n].k*c.dt) + Rint*c.dt2 ;
#ifdef CHECKS
        if(isnan(M[n1+m].R)) {
            printf("R[%d] = %lf in CR_step at time %d!!!\n",n1+m,M[n1+m].R,n);
            exit(EXIT_FAILURE) ;
        }
#endif
        //for(l=0;l<m;l++) C2int += (v[n-l].G+M[n0+l].MC)*M[tri(m,l)].R ;
        for(l=0;l<m;l++) C2int += M[n0+l].MC*M[tri(m,l)].R ;
        M[n1+m].C = M[n0+m].C*(1.-v[n].k*c.dt) + (C1int+C2int)*c.dt2 ;
#ifdef CHECKS
        if(isnan(M[n1+m].C)) {
            printf("C[%d] = %lf in CR_step at time %d!!!\n",n1+m,M[n1+m].C,n);
            exit(EXIT_FAILURE) ;
        }
#endif
    }
    
    M[nn].R = 0.5 ;
    double C1int=0., C2int=0. ;
    for(l=1;l<=n;l++) C1int += M[n0+l].MR*M[n1+l].C ;
    //for(l=0;l<=n;l++) C2int += (v[n-l].G+M[n0+l].MC)*M[n1+l].R ;
    for(l=0;l<=n;l++) C2int += M[n0+l].MC*M[n1+l].R ;  // jamming, no colored noise
    //printf("C1int = %lf\tC2int = %lf\tdt2 = %lf\n",C1int,C2int,c.dt2);
    M[nn].C = M[nn-1].C*(1.-v[n].k*c.dt) + (C1int+C2int)*c.dt2 ;
    //printf("C[%d] = %lf\n",nn,M[nn].C);
#ifdef CHECKS
    if(isnan(M[nn].C)) {
        printf("R[%d] = %lf in CR_step at time %d!!!\n",nn,M[nn].C,n);
        exit(EXIT_FAILURE) ;
    }
#endif
    
    v[n+1].chiR = v[n].chiR + (M[n0].R+M[n1].R)*c.dt/2. ;
}


// Kernel evaluation through integral over h0

void h0int(dyn_t *D, hint_t hi, coeffs_t *c, int n_out, int n, int nh0) {
    int nh, m, n1 = tri(n+1,0), n0 = tri(n,0) ; //n0 = tri(n,0), , nn=tri(n+1,n+1) ;
    int measure=n_out-n+1 ;
    c->sM0 = 1./sqrt(D->M[tri(n,n)].Vn);
    increments_t *i1 = (increments_t*) calloc(1,sizeof(increments_t)) ;
    i1->MR1 = (double*) calloc(n+2,sizeof(double)) ;
    i1->MC1 = (double*) calloc(n+2,sizeof(double)) ;
    
    FILE *fk ; fk = fopen("k_h0.dat","a");
    
    for(nh=nh0;nh<hi.NH;nh++)
    {
        i1->k1=0. ; i1->E1=0. ; i1->P1=0. ;
        for(m=0;m<n+2;m++) { i1->MR1[m]=0. ; i1->MC1[m]=0. ; }
        
        traj_step(D,i1,hi,*c,D->hw[nh].h0,measure,n,nh) ;
        //printf("nh = %04d\twh0*k1 = %21.15e\n",nh,D->hw[nh].wh0*i1.k1);
#ifdef CHECKS
        if(isnan(i1->k1)) {
            printf("k1 = %lf in h0int at time %d and nh=%d!!!\n",i1->k1,n,nh);
            exit(EXIT_FAILURE) ;
        }
        for(m=0;m<n+2;m++) {
            if(isnan(i1->MR1[m])) {
                printf("MR1[%d] = %lf in h0int at time %d and nh=%d!!!\n",m,i1->MR1[m],n,nh);
                exit(EXIT_FAILURE) ;
            }
            if(isnan(i1->MC1[m])) {
                printf("MC1[%d] = %lf in h0int at time %d and nh=%d!!!\n",m,i1->MC1[m],n,nh);
                exit(EXIT_FAILURE) ;
            }
        }
#endif
        if(measure) {
            D->v[n+1].k += D->hw[nh].wh0*i1->k1 ;
            fprintf(fk,"%lf\t%21.15e\n",D->hw[nh].h0,D->hw[nh].wh0*i1->k1) ;
            for(m=0;m<n+2;m++) D->M[n1+m].MR += D->hw[nh].wh0*i1->MR1[m] ;
            for(m=0;m<n+2;m++) D->M[n1+m].MC += D->hw[nh].wh0*i1->MC1[m] ;
            D->v[n+1].E += D->hw[nh].wh0*i1->E1 ;
            D->v[n+1].P += D->hw[nh].wh0*i1->P1 ;
        }
    }
    if(measure) {
        double c1 ;
        nh-- ;
        c1 = D->hw[nh].wh0*i1->k1/D->v[n+1].k ;
        for(m=0;m<n+2;m++) c1 = fmax(c1,D->hw[nh].wh0*i1->MR1[m]/D->M[n1+m].MR) ;
        for(m=0;m<n+2;m++) c1 = fmax(c1,D->hw[nh].wh0*i1->MC1[m]/D->M[n1+m].MC) ;
        D->v[n+1].err = c1 ;
        for(int m=0;m<=n;m++) {
            D->M[n1+m].Gtest *= hi.normNR ;
            D->M[n1+m].MCtest *= hi.normNR ;
        }
    }
    free(i1->MR1) ; free(i1->MC1) ; free(i1) ;
    fprintf(fk,"\n\n");
    fclose(fk) ;
    
    D->v[n+1].chiMR = D->v[n].chiMR + (D->M[n0].MR+D->M[n1].MR)*c->dt/2. ;
    
}


// Stochastic evolution of microscopic observable at given h0

void traj_step(dyn_t *D, increments_t *I1, hint_t hi, coeffs_t c, double h0, int measure, int n, int nh) {
    int np ;
    int n0=tri(n,0), n1=tri(n+1,0), nn=tri(n+1,n+1), l, m ;
    for(np=0;np<hi.NP;np++)
    {
        int nr = hi.NP*nh + np ;

        // Correlated noise generation
        
        D->X[n][nr].mu=0. ;
        for(m=c.m0;m<n;m++) D->X[n][nr].mu += D->M[n0+m].Vn*D->X[m][nr].eta ;
        D->X[n][nr].mu /= (-D->M[n0+n].Vn) ;
        D->X[n][nr].eta = D->X[n][nr].mu + c.sM0*GAUSSGEN ;
        if(measure==2) {
            for(m=0;m<=n;m++) {
                D->M[n0+m].MCtest += D->X[n][nr].eta*D->X[m][nr].eta ;
            }
        }
        
        // Evolution of y
        double Yint=0. ;
        for(l=1;l<=n;l++) Yint += D->M[n0+l].MR*D->X[l][nr].y ;
        //D->X[n+1][nr].y = D->X[n][nr].y*(1-D->v[n].k*c.dt) + Yint*c.dt2 +(D->X[n][nr].f+D->X[n][nr].eta-D->X[n][nr].v[1])*c.dt ;
        D->X[n+1][nr].y = D->X[n][nr].y*(1-D->v[n].k*c.dt) + Yint*c.dt2 +(D->X[n][nr].eta-D->X[n][nr].v[1])*c.dt ; // Jamming: no active forces
        
        // Evolution of h
        double h1 = h0 + D->X[n+1][nr].y + D->M[nn].C ;
        D->X[n+1][nr].h = h1 ;
        for(m=0;m<3;m++) D->X[n+1][nr].v[m] = c.eps*V(h1,m) ;
        
        // Evolution of H
        double Hint=0. ;
        for(m=0;m<=n;m++) {
            for(l=m+1;l<=n;l++) Hint += D->M[n0+l].MR*D->H[tri(l,m)][nr] ;
            D->H[n1+m][nr] = D->H[n0+m][nr]*(1.-(D->v[n].k+D->X[n][nr].v[2])*c.dt) + Hint*c.dt2 ;
        }
        D->H[nn][nr] = D->X[n+1][nr].v[2] ;
        
        // Active process
        /*
        D->X[n+1][nr].f = c.fc3*(c.fc1*D->X[n][nr].f+c.fc2*GAUSSGEN) ;
        if(measure==2) {
            for(m=0;m<=n+1;m++) {
                D->M[n1+m].Gtest += D->X[n+1][nr].f*D->X[m][nr].f ;
            }
        }
       */
        
        // Kernel update
        if(measure) {
             I1->k1 += D->X[n+1][nr].v[2] + D->X[n+1][nr].v[1] ;
#ifdef CHECKS
            if(isnan(I1->k1)) {
                printf("I1->k1 = %lf in traj_step at time %d, nh=%d and np=%d!!!\n",I1->k1,n,nh,np);
                printf("y[%d][%d] = %lf\n",n,nr,D->X[n][nr].y);
                printf("h[%d][%d] = %lf\n",n,nr,D->X[n][nr].h);
                printf("eta[%d][%d] = %lf\n",n,nr,D->X[n][nr].eta);
                printf("f[%d][%d] = %lf\n",n,nr,D->X[n][nr].f);
                for(m=0;m<3;m++) printf("v[%d][%d][%d] = %lf\t",m,n,nr,D->X[n][nr].v[m]);
                printf("\n");
                printf("sM0 = %lf\n",c.sM0);
                printf("Vn[%d] = %lf\n",tri(n,n),D->M[tri(n,n)].Vn);
                if(n<=10) {
                    double **A = (double**) calloc(n+1,sizeof(double*));
                    for(m=0;m<n+1;m++) {
                        A[m] = (double*) calloc(n+2,sizeof(double)) ;
                        for(l=0;l<n+1;l++) A[m][l] = D->M[tri(m,l)].MC ;
                    }
                    printf("det MC = %21.15e\n",det(A,n+1)) ;
                    for(m=0;m<n+2;m++) free(A[m]) ;
                    free(A) ;
                }
                
                FILE *f ; f = fopen("MC.dat","w") ;
                int nmax = tri(n,n) ;
                for(m=0;m<=nmax;m++) fprintf(f,"%21.15e\n",D->M[m].MC) ;
                fclose(f) ;
                
                //phys_t par ; par.dt = c.dt ;
                //print_output(D->v,D->M,par) ;
                
                exit(EXIT_FAILURE) ;
            }
#endif
            for(m=0;m<n+2;m++) I1->MR1[m] += D->X[n+1][nr].v[2]*D->H[n1+m][nr] ;
            for(m=0;m<n+2;m++) I1->MC1[m] += D->X[n+1][nr].v[1]*D->X[m][nr].v[1] ;
            I1->E1 += D->X[n+1][nr].v[0] ;
            I1->P1 -= D->X[n+1][nr].v[1] ;
        }
    }
}



// MATRIX INVERSION FUNCTIONS //


// Matrix inversion with blocks scheme.

void block_matrix_inversion(mat_t *M, int n, int m0) {
    
    int m, l, n1=tri(n+1,0), nn=tri(n+1,n+1) ;
    double *vM = (double*) calloc(n+1,sizeof(double)) ;
    for(m=m0;m<=n;m++) for(l=0;l<=n;l++) vM[m] += M[tri(m,l)].Minv*M[n1+l].MC ;
    
    double sigma2M = 0. ;
    for(m=m0;m<=n;m++) sigma2M += M[n1+m].MC*vM[m] ;
    sigma2M = 1./(M[nn].MC-sigma2M) ;
    
    for(m=0;m<=n;m++)
    {
        for(l=0;l<=m;l++) M[tri(m,l)].Minv += sigma2M*vM[m]*vM[l] ;
        M[n1+m].Minv = -sigma2M*vM[m] ;
        M[n1+m].Vn = M[n1+m].Minv ;
    }
    M[nn].Minv = sigma2M ;
    M[nn].Vn = M[nn].Minv ;
    gsl_matrix *MM = gsl_matrix_calloc(n+2,n+2) ;
    gsl_matrix *MI = gsl_matrix_calloc(n+2,n+2) ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            gsl_matrix_set(MM,m,l,M[tri(m,l)].MC) ;
            gsl_matrix_set(MI,m,l,M[tri(m,l)].Minv) ;
        }
    }   
    FILE *fe ; fe = fopen("err.dat","a") ;
    fprintf(fe,"%03d\t%21.15e\n",n+2,check_id_gsl(MI,MM,n+2)) ;
    fclose(fe) ;
    gsl_matrix_free(MM) ;
    gsl_matrix_free(MI) ;
    free(vM) ;
}


// Compute the average of the last n_err (or less) elements of Minv * M ; the result is the deviation from the correct inverse matrix

double check_id(mat_t *M, int n, int n_err) {
    
    int m, l, m0=n+2-n_err ;
    m0 = m0*(m0>0) ;
    double sum_err=0., norm = (double) (n+1-m0) ;
    for(m=m0+1;m<n+2;m++) {
        double x=0. ;
        for(l=0;l<n+2;l++) x += M[tri(m,l)].Minv*M[tri(l,0)].MC ;
        sum_err += fabs(x) ;
    }
    sum_err /= norm ;
    return sum_err ;
}


// Invert the matrix MC at step n using LU decomposition and inversion in GSL library, and replace Minv with the new inverse


void stabilize_cholesky(mat_t *M, int n) {
 
    int m, l, n1=tri(n+1,0) ;
    
#ifdef CHECK_MATRICES
    FILE *f0 ; char file0[FLEN]={'\0'};
    sprintf(file0,"M_%d.dat",n) ; f0 = fopen(file0,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f0,"%21.15e\t",M[tri(m,l)].MC) ;
        }
        fprintf(f0,"\n");
    }
    fclose(f0) ;
    FILE *f1 ; char file1[FLEN]={'\0'};
    sprintf(file1,"MI_%d.dat",n) ; f1 = fopen(file1,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f1,"%21.15e\t",M[tri(m,l)].Minv) ;
        }
        fprintf(f1,"\n");
    }
    fclose(f1) ;
#endif
    
    gsl_matrix *MM = gsl_matrix_calloc(n+2,n+2) ;
    gsl_matrix *MMI = gsl_matrix_calloc(n+2,n+2) ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            gsl_matrix_set(MM,m,l,M[tri(m,l)].MC) ;
            gsl_matrix_set(MMI,m,l,M[tri(m,l)].Minv) ;
        }
    }
    
#ifdef CHECK_MATRICES
    gsl_matrix *I = gsl_matrix_calloc(n+2,n+2) ;
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,MMI,MM,0.,I);
    FILE *f1I ; char file1I[FLEN]={'\0'};
    sprintf(file1I,"I_%d.dat",n) ; f1I = fopen(file1I,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f1I,"%21.15e\t",gsl_matrix_get(I,m,l)) ;
        }
        fprintf(f1I,"\n");
    }
    fclose(f1I) ;
#endif
    
#ifdef STABILIZATION
    gsl_matrix_set_zero(MMI) ;
    gsl_matrix *MCH = gsl_matrix_calloc(n+2,n+2) ;
    
    gsl_matrix_memcpy(MCH,MM) ;
    gsl_linalg_cholesky_decomp1(MCH) ;
    gsl_linalg_cholesky_invert(MCH) ;
    gsl_matrix_memcpy(MMI,MCH) ;
    
    gsl_matrix *R = gsl_matrix_calloc(n+2,n+2) ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<=m;l++) {
            gsl_matrix_set(R,m,l,M[tri(m,l)].Minv/gsl_matrix_get(MMI,m,l)-1.) ;
            M[tri(m,l)].Minv = gsl_matrix_get(MMI,m,l) ;
        }
    }
    for(m=0;m<n+2;m++) M[n1+m].Vn = M[n1+m].Minv ;
    gsl_matrix_free(MCH) ;
    
    printf("Stabilization done at n=%d\t",n) ;
    printf("err(%d) = %21.15e\n",n,check_id(M,n,5));
    
#ifdef CHECK_MATRICES
    
    FILE *f2 ; char file2[FLEN]={'\0'};
    sprintf(file2,"MILU_%d.dat",n) ; f2 = fopen(file2,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f2,"%21.15e\t",M[tri(m,l)].Minv) ;
        }
        fprintf(f2,"\n");
    }
    fclose(f2) ;
    FILE *fR ; char fileR[FLEN]={'\0'};
    sprintf(fileR,"R_%d.dat",n) ; fR = fopen(fileR,"w") ;
    gsl_matrix_fprintf(fR,R,"%e") ;
    fclose(fR) ;
    gsl_matrix_set_zero(I) ;
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,MMI,MM,0.,I);
    FILE *f2I ; char file2I[FLEN]={'\0'};
    sprintf(file2I,"ILU_%d.dat",n) ; f2I = fopen(file2I,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f2I,"%21.15e\t",gsl_matrix_get(I,m,l)) ;
        }
        fprintf(f2I,"\n");
    }
    fclose(f2I) ;
    
    gsl_matrix_free(I) ;
#endif
    gsl_matrix_free(R) ;
#endif
    
    gsl_matrix_free(MMI) ;
    gsl_matrix_free(MM) ;
    
}


void Cholesky_inversion(mat_t *M, int n) {
    
    int m, l, n1=tri(n+1,0) ;
    gsl_matrix *MM = gsl_matrix_calloc(n+2,n+2) ;
    gsl_matrix *MC = gsl_matrix_calloc(n+2,n+2) ;
    for(m=0;m<n+2;m++) for(l=0;l<n+2;l++) gsl_matrix_set(MC,m,l,M[tri(m,l)].MC) ;
    gsl_matrix_memcpy(MM,MC) ;
    gsl_linalg_cholesky_decomp1(MM) ;
    gsl_linalg_cholesky_invert(MM) ;
    for(m=0;m<n+2;m++) for(l=0;l<=m;l++) M[tri(m,l)].Minv = gsl_matrix_get(MM,m,l) ;
    for(m=0;m<n+2;m++) M[n1+m].Vn = M[n1+m].Minv ;
    FILE *fe ; fe = fopen("err.dat","a") ;
    fprintf(fe,"%03d\t%21.15e\n",n+2,check_id_gsl(MM,MC,n+2)) ;
    fclose(fe) ;
    gsl_matrix_free(MM) ;
    gsl_matrix_free(MC) ;
}

void LU_inversion(mat_t *M, int n) {
    int sgn, m, l, n1=tri(n+1,0) ;
    gsl_matrix *MM = gsl_matrix_calloc(n+2,n+2) ;
    gsl_matrix *MLU = gsl_matrix_calloc(n+2,n+2) ;
    gsl_matrix *MI = gsl_matrix_calloc(n+2,n+2) ;
    gsl_permutation *p = gsl_permutation_calloc(n+2) ;
    for(m=0;m<n+2;m++) for(l=0;l<n+2;l++) gsl_matrix_set(MM,m,l,M[tri(m,l)].MC) ;
    gsl_matrix_memcpy(MLU,MM) ;
    gsl_linalg_LU_decomp(MLU,p,&sgn) ;
    gsl_linalg_LU_invert(MLU,p,MI) ;
    for(m=0;m<n+2;m++) for(l=0;l<=m;l++) M[tri(m,l)].Minv = gsl_matrix_get(MI,m,l) ;
    for(m=0;m<n+2;m++) M[n1+m].Vn = M[n1+m].Minv ;
    FILE *fe ; fe = fopen("err.dat","a") ;
    fprintf(fe,"%03d\t%21.15e\n",n+2,check_id_gsl(MM,MI,n+2)) ;
    fclose(fe) ;
    gsl_matrix_free(MM) ;
    gsl_matrix_free(MLU) ;
    gsl_matrix_free(MI) ;
    gsl_permutation_free(p) ;
}

double check_id_gsl(gsl_matrix *A, gsl_matrix *B, int n) {
    int i, j;
    double err=0. ;
    gsl_matrix *I1 = gsl_matrix_calloc(n,n) ;
    gsl_matrix *I2 = gsl_matrix_calloc(n,n) ;
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,A,B,0.,I1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,B,A,0.,I2);
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            //printf("err(%d,%d) = %21.15e\n",i,j,gsl_matrix_get(I,i,j)) ;
            err = fmax(err,fabs(gsl_matrix_get(I1,i,j)-(double)((i==j)))) ;
            err = fmax(err,fabs(gsl_matrix_get(I2,i,j)-(double)((i==j)))) ;
        }
    }
    gsl_matrix_free(I1) ;
    gsl_matrix_free(I2) ;
    return err ;
}



void stabilize(mat_t *M, int n) {
    
    int sgn, m, l, n1=tri(n+1,0) ;
    
#ifdef CHECK_MATRICES
    FILE *f0 ; char file0[FLEN]={'\0'};
    sprintf(file0,"M_%d.dat",n) ; f0 = fopen(file0,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f0,"%21.15e\t",M[tri(m,l)].MC) ;
        }
        fprintf(f0,"\n");
    }
    fclose(f0) ;
    FILE *f1 ; char file1[FLEN]={'\0'};
    sprintf(file1,"MI_%d.dat",n) ; f1 = fopen(file1,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f1,"%21.15e\t",M[tri(m,l)].Minv) ;
        }
        fprintf(f1,"\n");
    }
    fclose(f1) ;
#endif
    
    gsl_matrix *MM = gsl_matrix_calloc(n+2,n+2) ;
    gsl_matrix *MMI = gsl_matrix_calloc(n+2,n+2) ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            gsl_matrix_set(MM,m,l,M[tri(m,l)].MC) ;
            gsl_matrix_set(MMI,m,l,M[tri(m,l)].Minv) ;
        }
    }
    
#ifdef CHECK_MATRICES
    gsl_matrix *I = gsl_matrix_calloc(n+2,n+2) ;
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,MMI,MM,0.,I);
    FILE *f1I ; char file1I[FLEN]={'\0'};
    sprintf(file1I,"I_%d.dat",n) ; f1I = fopen(file1I,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f1I,"%21.15e\t",gsl_matrix_get(I,m,l)) ;
        }
        fprintf(f1I,"\n");
    }
    fclose(f1I) ;
#endif
    
#ifdef STABILIZATION
    gsl_matrix_set_zero(MMI) ;
    gsl_matrix *MLU = gsl_matrix_calloc(n+2,n+2) ;
    gsl_permutation *p = gsl_permutation_calloc(n+2) ;
    
    gsl_matrix_memcpy(MLU,MM) ;
    gsl_linalg_LU_decomp(MLU,p,&sgn) ;
    gsl_linalg_LU_invert(MLU,p,MMI) ;
    
    gsl_matrix *R = gsl_matrix_calloc(n+2,n+2) ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<=m;l++) {
            gsl_matrix_set(R,m,l,M[tri(m,l)].Minv/gsl_matrix_get(MMI,m,l)-1.) ;
            M[tri(m,l)].Minv = gsl_matrix_get(MMI,m,l) ;
        }
    }
    for(m=0;m<n+2;m++) M[n1+m].Vn = M[n1+m].Minv ;
    gsl_matrix_free(MLU) ;
    gsl_permutation_free(p) ;
    
    printf("Stabilization done at n=%d\t",n) ;
    printf("err(%d) = %21.15e\n",n,check_id(M,n,5));
    
#ifdef CHECK_MATRICES
    
    FILE *f2 ; char file2[FLEN]={'\0'};
    sprintf(file2,"MILU_%d.dat",n) ; f2 = fopen(file2,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f2,"%21.15e\t",M[tri(m,l)].Minv) ;
        }
        fprintf(f2,"\n");
    }
    fclose(f2) ;
    FILE *fR ; char fileR[FLEN]={'\0'};
    sprintf(fileR,"R_%d.dat",n) ; fR = fopen(fileR,"w") ;
    gsl_matrix_fprintf(fR,R,"%e") ;
    fclose(fR) ;
    gsl_matrix_set_zero(I) ;
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,MMI,MM,0.,I);
    FILE *f2I ; char file2I[FLEN]={'\0'};
    sprintf(file2I,"ILU_%d.dat",n) ; f2I = fopen(file2I,"w") ;
    for(m=0;m<n+2;m++) {
        for(l=0;l<n+2;l++) {
            fprintf(f2I,"%21.15e\t",gsl_matrix_get(I,m,l)) ;
        }
        fprintf(f2I,"\n");
    }
    fclose(f2I) ;
    
    gsl_matrix_free(I) ;
#endif
    gsl_matrix_free(R) ;
#endif
    
    gsl_matrix_free(MMI) ;
    gsl_matrix_free(MM) ;
}






// Reallocation of new trajectories when NH <- NH+1

void resize(dyn_t *D, hint_t *hi, int N) {
    
    hi->NH += 1 ; int NR = hi->NH*hi->NP, n, N1=(int)(N*(N+1)/2) ;
    hweights_t *tmp_hw = realloc(D->hw,hi->NH*sizeof(hweights_t)) ;
    if(tmp_hw!=NULL) D->hw = tmp_hw ;
    else { printf("Reallocation of hweights failed!!\n\n"); exit(EXIT_FAILURE) ; }
    for(n=0;n<N1;n++)
    {
        double *tmp_H = realloc(D->H[n],NR*sizeof(double));
        if(tmp_H!=NULL) D->H[n] = tmp_H ;
        else { printf("Reallocation of H[%d] failed!!\n\n",n); exit(EXIT_FAILURE) ; }
    }
    for(n=0;n<N;n++)
    {
        traj_t *tmp_X = realloc(D->X[n],NR*sizeof(traj_t));
        if(tmp_X!=NULL) D->X[n] = tmp_X ;
        else { printf("Reallocation of X[%d] failed!!\n\n",n); exit(EXIT_FAILURE) ; }
    }
}


// Normalization of numerically measured correlations

void norm_test(dyn_t *D, hint_t hi, int N) {
    double NRd = ((double)(hi.NH*hi.NP)) ;
    for(int n=0;n<N;n++) {
        for(int m=0;m<=n;m++) {
            D->M[tri(n,m)].Gtest /= NRd ;
            D->M[tri(n,m)].MCtest /= NRd ;
        }
    }
}





// Standard for output file names

void filename(char *s, char *v, phys_t par) {
//sprintf(s,"%s_%.2f_%.0e_%.0e_%.0e_%.0e.dat",v,par.phi,par.eps,par.Teff,par.b0,par.taup);
sprintf(s,"%s_%.2f_%.0e_%.0e.dat",v,par.phi,par.eps,par.b0);
}



// Print of all the outputs files

void print_output(vec_t *v, mat_t *M, phys_t par) {
    
    FILE *fC, *fR, *fk, *fMC, *fMR, *fMCt /*, *fG, *fGt*/ ;
    FILE *fD, *fdR, *fdMC, *fdMR, *fE, *fP, *fz, *ferr ;
    char fileC[FLEN]={'\0'}; filename(fileC,"C",par) ;
    char fileR[FLEN]={'\0'}; filename(fileR,"R",par) ;
    char fileD[FLEN]={'\0'}; filename(fileD,"D",par) ;
    char filek[FLEN]={'\0'}; filename(filek,"k",par) ;
    char fileMC[FLEN]={'\0'}; filename(fileMC,"MC",par) ;
    char fileMR[FLEN]={'\0'}; filename(fileMR,"MR",par) ;
    //char fileG[FLEN]={'\0'}; filename(fileG,"G",par) ;
    char fileMCt[FLEN]={'\0'}; filename(fileMCt,"MCt",par) ;
    //char fileGt[FLEN]={'\0'}; filename(fileGt,"Gt",par) ;
    char filedR[FLEN]={'\0'}; filename(filedR,"dR",par) ;
    char filedMC[FLEN]={'\0'}; filename(filedMC,"dMC",par) ;
    char filedMR[FLEN]={'\0'}; filename(filedMR,"dMR",par) ;
    char fileE[FLEN]={'\0'}; filename(fileE,"E",par) ;
    char fileP[FLEN]={'\0'}; filename(fileP,"P",par) ;
    char filez[FLEN]={'\0'}; filename(filez,"z",par) ;
    char filerr[FLEN]={'\0'}; filename(filerr,"err",par) ;
    
    fC = fopen(fileC,"w") ;
    fR = fopen(fileR,"w") ;
    fD = fopen(fileD,"w") ;
    fk = fopen(filek,"w") ;
    fMC = fopen(fileMC,"w") ;
    fMR = fopen(fileMR,"w") ;
    //fG = fopen(fileG,"w") ;
    fMCt = fopen(fileMCt,"w") ;
    //fGt = fopen(fileGt,"w") ;
    fdR = fopen(filedR,"w") ;
    fdMC = fopen(filedMC,"w") ;
    fdMR = fopen(filedMR,"w") ;
    fE = fopen(fileE,"w") ;
    fP = fopen(fileP,"w") ;
    fz = fopen(filez,"w") ;
    ferr = fopen(filerr,"w") ;
    
    int n,m ;
    for(n=0;n<par.N;n++) {
        double t = n*par.dt ;
        int nn = tri(n,n) ;
        fprintf(fk,"%lf\t%21.15e\t%21.15e\t%21.15e\n",t,v[n].k,v[n].chiR,v[n].chiMR);
        //fprintf(fG,"%lf\t%21.15e\n",t,v[n].G);
        fprintf(fE,"%lf\t%21.15e\n",t,v[n].E);
        fprintf(fP,"%lf\t%21.15e\n",t,v[n].P);
	fprintf(fz,"%lf\t%21.15e\n",t,v[n].k+v[n].P);
        fprintf(ferr,"%lf\t%21.15e\n",t,v[n].err);
        fprintf(fdR,"%lf\t%21.15e\n",t,M[nn].R);
        fprintf(fdMC,"%le\t%21.15e\n",t,M[nn].MC);
        fprintf(fdMR,"%lf\t%21.15e\n",t,M[nn].MR);
        fprintf(fC,"%lf\t",t);
        fprintf(fR,"%lf\t",t);
        fprintf(fD,"%lf\t",t);
        fprintf(fMC,"%lf\t",t);
        fprintf(fMR,"%lf\t",t);
        fprintf(fMCt,"%lf\t",t);
        //fprintf(fGt,"%lf\t",t);
        
        for(m=0;m<=n;m++) {
            int nm = tri(n,m) ;
            fprintf(fC,"%lf\t",M[nm].C);
            fprintf(fR,"%lf\t",M[nm].R);
            fprintf(fD,"%lf\t",M[nn].C+M[tri(m,m)].C-2.*M[nm].C);
            fprintf(fMC,"%lf\t",M[nm].MC);
            fprintf(fMR,"%lf\t",M[nm].MR);
            fprintf(fMCt,"%lf\t",M[nm].MCtest);
            //fprintf(fGt,"%lf\t",M[nm].Gtest);
        }
        for(m=n+1;m<par.N;m++) {
            int nm = tri(n,m) ;
            fprintf(fC,"%lf\t",M[nm].C);
            fprintf(fR,"%lf\t",0.);
            fprintf(fD,"%lf\t",M[nn].C+M[tri(m,m)].C-2.*M[nm].C);
            fprintf(fMC,"%lf\t",M[nm].MC);
            fprintf(fMR,"%lf\t",0.);
            fprintf(fMCt,"%lf\t",M[nm].MCtest);
            //fprintf(fGt,"%lf\t",M[nm].Gtest);
        }
        fprintf(fC,"\n");
        fprintf(fR,"\n");
        fprintf(fD,"\n");
        fprintf(fMC,"\n");
        fprintf(fMR,"\n");
        fprintf(fMCt,"\n");
        //fprintf(fGt,"\n");
        
    }
    fclose(fk) ; //fclose(fG) ;
    fclose(fdR) ; fclose(fdMC) ; fclose(fdMR) ;
    fclose(fC) ; fclose(fR) ; fclose(fD) ;
    fclose(fMC) ; fclose(fMR) ;
    fclose(fMCt) ; //fclose(fGt) ;
    fclose(fE) ; fclose(fP) ; fclose(fz) ;
    fclose(ferr) ;

}





void print_info(phys_t par, hint_t hi, ini_t I0) {
 
    FILE *fp ;
    fp = fopen("INFO","w") ;
    
    fprintf(fp,"Physical parameters :\n\n");
    fprintf(fp," phi = %.2e\n Teff = %.2e\n eps = %.2e\n taup = %.2e\n\n",par.phi,par.Teff,par.eps,par.taup) ;
    fprintf(fp," Time step = %.2e\n Steps = %d\n Total time = %.2e\n\n\n",par.dt,par.N,par.N*par.dt);
    
    fprintf(fp,"Integration over h0 :\n\n") ;
    fprintf(fp," hmin = %lf\n hmax = %lf\n dh = %lf\n NH = %d\n z = %lf\n\n",hi.hmin,hi.hmax,hi.dh0,hi.NH,hi.z) ;
    fprintf(fp," Threshold for the relative error = %.2e\n\n",hi.threshold) ;
    fprintf(fp," Number of realizations for each h0 : NP = %d\n\n",hi.NP) ;
//fprintf(fp," Random generator seed : SEED = %d\n\n\n",SEED) ;
    
    fprintf(fp,"Initial conditions :\n\n") ;
    fprintf(fp," k0 = %.6e\n",I0.k0) ;
    fprintf(fp," MR0 = %.6e\n",I0.MR0) ;
    fprintf(fp," MC0 = %.6e\n",I0.MC0) ;
    fprintf(fp," E0 = %.6e\n",I0.E0) ;
    fprintf(fp," P0 = %.6e\n",I0.P0) ;
    fprintf(fp," G0 = %.6e\n",I0.G0) ;
    
    fclose(fp) ;
}


void free_dyn(dyn_t *D, int N, int N1) {
    int n ;
    for(n=0;n<N;n++) free(D->X[n]) ;
    for(n=0;n<N1;n++) free(D->H[n]) ;
    free(D->X) ; free(D->H) ; free(D->v) ; free(D->M) ; free(D->hw) ;
    free(D) ;
}

void print_coeffs(coeffs_t c) {
 
    FILE *fp ;
    fp = fopen("COEFFS","w") ;

    fprintf(fp,"m0 = %d\n",c.m0) ;
    
    fclose(fp) ;
}



 // Recursive determinant computator

 double det(double **M, int n) {
     
     if(n==1) return M[0][0] ;
     else {
         int j, k, l, m, s=1 ;
         double D=0. ;
         for(m=0;m<n;m++) {
             double** Mm = (double**) calloc(n-1,sizeof(double*)) ;
             for(l=0;l<n-1;l++) Mm[l] = (double*) calloc(n-1,sizeof(double)) ;
             j=0 ;
             for(l=0;l<n;l++) {
                 if(l==m) continue ;
                 else {
                     for(k=0;k<n-1;k++) Mm[j][k] = M[l][k+1] ;
                     j++ ;
                 }
             }
             D += s*M[m][0]*det(Mm,n-1) ;
             for(l=0;l<n-1;l++) free(Mm[l]) ;
             free(Mm) ;
             s *= -1 ;
         }
         return D ;
     }
 }


 // Matrix inversion with cofactors

 void matrix_inversion(mat_t *M, int n) {
     
     int j,k,l,m, s=1, sr=1 ;
 #ifdef CHECK_MATRICES
     printf("\n\t\t MC, n = %d\n\n",n);
     for(m=0;m<n;m++) {
         for(l=0;l<=m;l++) printf("%21.15e\t",M[tri(m,l)].MC);
         printf("\n");
     }
 #endif
     
     
     double **A = (double**) calloc(n,sizeof(double*)) ;
     double **B = (double**) calloc(n,sizeof(double*)) ;
     for(m=0;m<n;m++) {
         A[m] = (double*) calloc(n,sizeof(double)) ;
         B[m] = (double*) calloc(n,sizeof(double)) ;
         for(l=0;l<n;l++) A[m][l] = M[tri(m,l)].MC ;
     }

     /*
 #ifdef CHECKS
     printf("\n\t\t A = MC ?, n = %d\n\n",n);
     for(m=0;m<n;m++) {
         for(l=0;l<n;l++) printf("%21.15e\t",A[m][l]);
         printf("\n");
     }
 #endif
      */
     
     double D = det(A,n) ;
 //#ifdef CHECKS
     printf("\ndet A(%d) = %21.15e\n",n,D) ;
 //#endif
     
     for(m=0;m<n;m++) {
         s=sr ;
         for(l=0;l<n;l++) {
             double **Mm = (double**) calloc(n-1,sizeof(double*)) ;
             for(j=0;j<n-1;j++) Mm[j] = (double*) calloc(n-1,sizeof(double)) ;
             int j1=0;
             for(j=0;j<n;j++) {
                 int k1=0 ;
                 if(j==m) continue ;
                 else {
                     for(k=0;k<n;k++) {
                         if(k==l) continue ;
                         else {
                             Mm[j1][k1] = A[j][k] ;
                             k1++ ;
                         }
                     }
                     j1++ ;
                 }
             }
             //printf("det Mm = %lf \t",det(Mm,n-1)) ;
             B[l][m] = s*det(Mm,n-1)/D ;
             M[tri(l,m)].Minv = B[l][m] ;
             for(j=0;j<n-1;j++) free(Mm[j]) ;
             free(Mm) ;
             s *= -1 ;
         }
         sr *= -1 ;
     }
     for(m=0;m<n;m++) M[tri(n-1,m)].Vn = M[tri(n-1,m)].Minv ;
     
 #ifdef CHECK_MATRICES
     /*
     printf("\n\t\t B, n = %d\n\n",n);
     for(m=0;m<n;m++) {
         for(l=0;l<n;l++) {
             printf("%21.15e\t",B[m][l]);
         }
         printf("\n");
     }
      */
     
     
     printf("\n\t\t Minv, n = %d\n\n",n);
     for(m=0;m<n;m++) {
         for(l=0;l<n;l++) {
             printf("%21.15e\t",M[tri(m,l)].Minv);
         }
         printf("\n");
     }
     
     printf("\n\t\t Vn, n = %d\n\n",n) ;
     for(m=0;m<n;m++) {
         for(l=0;l<=m;l++) printf("%21.15e\t",M[tri(m,l)].Vn);
         printf("\n");
     }
     printf("\n\n") ;
     printf("\n M * Minv = Id ?\n\n") ;
     for(m=0;m<n;m++) {
         for(l=0;l<n;l++) {
             double x = 0. ;
             for(int j=0;j<n;j++) x += M[tri(m,j)].Minv*M[tri(j,l)].MC ;
             printf("%21.15e\t",x) ;
         }
         printf("\n\n") ;
     }
 #endif
     
 }


/*
void print_output(vec_t*v, mat_t *M, phys_t par) {
 
    
    int U=7, S=5, R=2, u,s,r ;
    char **su = (char**) malloc(U*sizeof(char*))   ;
    //char **ss, *    *sr ;
    su[0] = "E" ; su[1] = "P" ; su[2] = "k" ; su[3] = "G" ;
    su[4] = "dR" ; su[5] = "dMR" ; su[6] = "dMC" ;
    //ss[0] = "C" ; ss[1] = "D" ; ss[2] = "MC" ; ss[3] = "MCt" ; ss[4] = "Gt" ;
    //sr[0] = "R" ; sr[1] = "MR";
    
    //FILE *fu[U] , *fs[S] , *fr[R] ;
    //char *fileu[U], *files[S],  *filer[R] ;
    //FILE **fu = (FILE**) malloc(U*sizeof(FILE*)) ;
    char **fileu = (char**) malloc(U*sizeof(char*)) ;
    
    for(u=0;u<U;u++) {
        //filename(fileu[u],su[u],par) ;
        printf("%s\t%s\n",su[u],fileu[u]) ;
        //fu[u] = fopen(fileu[u],"w") ;
    }
    /
    for(s=0;s<S;s++) { filename(files[s],ss[s],par) ; fs[s] = fopen(files[s],"w") ; }
    for(r=0;r<R;r++) { filename(filer[r],sr[r],par) ; fr[r] = fopen(filer[r],"w") ; }
    
    int n, m ;
    
    
    for(n=0;n<par.N;n++) {
        double t = n*par.dt ; int nn = tri(n,n), n0 = tri(n,0) ;
        
        double *un = (double*) calloc(U,sizeof(double)) ;
        un[0] = v[n].E ; un[1] = v[n].P ; un[2] = v[n].k ; un[3] = v[n].G ;
        un[4] = M[nn].R ; un[5] = M[nn].MR ; un[6] = M[nn].MC ;
        
        for(u=0;u<U;u++) fprintf(fu[u],"%lf\t%21.15e\n",t,un[u]) ;
        free(un) ;
        
        for(s=0;s<S;s++) fprintf(fs[s],"%lf\t",t) ;
        for(r=0;r<R;r++) fprintf(fr[r],"%lf\t",t) ;
        
        for(m=0;m<par.N;m++) {
            int nm = tri(n,m) ;
            double *sn = (double*) calloc(S,sizeof(double)) ;
            sn[0] = M[nm].C ; sn[1] = M[nn].C+M[tri(m,m)].C-2.*M[nm].C ;
            sn[2] = M[nm].MC ; sn[3] = M[nm].MCtest ;
            sn[4] = M[nm].Gtest ;
            for(s=0;s<S;s++) fprintf(fs[s],"%lf\t",sn[s]) ;
            free(sn) ;
        }
        for(m=0;m<=n;m++) {
            double *rn = (double*) calloc(R,sizeof(double)) ;
            rn[0] = M[n0+m].R ; rn[1] = M[n0+m].MR ;
            for(r=0;r<R;r++) fprintf(fr[r],"%lf\t",rn[r]) ;
            free(rn) ;
        }
        for(s=0;s<S;s++) fprintf(fs[s],"\n") ;
        for(r=0;r<R;r++) fprintf(fr[r],"\n") ;
    }
    
    //for(u=0;u<U;u++) fclose(fu[u]) ;
    
    for(s=0;s<S;s++) fclose(fs[s]) ;
    for(r=0;r<R;r++) fclose(fr[r]) ;
    
}
*/


// codice per stampare i potenziali nelle traiettorie
/*
FILE *fv0 ; fv0 = fopen("v0.dat","w") ;
fprintf(fv0,"%lf",0.);
for(int nh=0;nh<hi0.NH;nh++) fprintf(fv0,"\t%21.15e",D->X[0][nh*hi0.NP].v[0]) ;
fprintf(fv0,"\n");
FILE *fv1 ; fv1 = fopen("v1.dat","w") ;
fprintf(fv1,"%lf",0.);
for(int nh=0;nh<hi0.NH;nh++) fprintf(fv1,"\t%21.15e",D->X[0][nh*hi0.NP].v[1]) ;
fprintf(fv1,"\n");
FILE *fv2 ; fv2 = fopen("v2.dat","w") ;
fprintf(fv2,"%lf",0.);
for(int nh=0;nh<hi0.NH;nh++) fprintf(fv2,"\t%21.15e",D->X[0][nh*hi0.NP].v[2]) ;
fprintf(fv2,"\n");


       fprintf(fv0,"%lf",(n+1)*p0.dt);
       for(int nh=0;nh<hi0.NH;nh++) fprintf(fv0,"\t%21.15e",D->X[n+1][nh*hi0.NP].v[0]) ;
       fprintf(fv0,"\n");
       fprintf(fv1,"%lf",(n+1)*p0.dt);
       for(int nh=0;nh<hi0.NH;nh++) fprintf(fv1,"\t%21.15e",D->X[n+1][nh*hi0.NP].v[1]) ;
       fprintf(fv1,"\n");
       fprintf(fv2,"%lf",(n+1)*p0.dt);
       for(int nh=0;nh<hi0.NH;nh++) fprintf(fv2,"\t%21.15e",D->X[n+1][nh*hi0.NP].v[2]) ;
       fprintf(fv2,"\n");

fclose(fh); fclose(fv0); fclose(fv1); fclose(fv2);
 
 */


// vecchio codice con l'inversione
 /*
        while(D->v[n+1].err>hi0.threshold) {
            resize(D,&hi0,N) ;                             // Reallocation of memory for D
            init_traj(D,p0,hi0,sqrt(i0.G0),hi0.NH-1) ;     // Initialization of new trajectories
            for(int n1=0;n1<=n;n1++) h0int(D,hi0,&c0,n,n1,hi0.NH-1) ; // Evolution from t=0 to t=n*dt ; it doesn't update measures at previous times
            printf("n = %d\tNH = %d\t err = %21.15e\n",n,hi0.NH,D->v[n+1].err);
        }
         */
    
        //matrix_inversion(D->M,n+2) ;
        
        // Inversion of MC matrix using blocks
        
        /*
        Cholesky_inversion(D->M,n) ;
        double err=check_id(D->M,n,5) ;
        printf("err(%d) = %21.15e\t\t",n,err);
        printf("s2= %21.15e\n",D->M[tri(n+1,n+1)].Minv) ;
        */
        
        //Cholesky_inversion(D->M,n) ;
        //LU_inversion(D->M,n) ;

/*
double err=check_id(D->M,n,5) ;
printf("err(%d) = %21.15e\n",n,err);
if(err>err_max) {
    stabilize_cholesky(D->M,n) ;
    err=check_id(D->M,n,5) ;
}
 */
