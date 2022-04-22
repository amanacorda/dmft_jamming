//
//  jamming_iterative.c
//  
//
//  Created by Alessandro Manacorda on 1/27/21.
//

#include "jamming_it.h"


int main(int argc, char **argv) {
    
    phys_t p0 ; hint_t hi0 ; ini_t i0 ;         // Declaration of physical parameters, integration constants, initial conditions
    
    init_par(argc-1,argv,&p0,&hi0,&i0) ;        // Initialization of parameters
    coeffs_t c0 ;                               // Declaration of numerical coefficients
    init_coeffs(p0,i0,&c0) ;                    // Initialization
    
    dyn_t *D = malloc(sizeof(dyn_t)) ;          // Declaration of dynamical variables
    init_dyn_it(D,p0,&hi0,i0) ;                    // Initialization
    increments_t *K = malloc(sizeof(increments_t)) ;
    double err[2]={0.}, CONV=1., CONV_TH=hi0.threshold;
    
    gsl_vector *eval = gsl_vector_alloc(N);
    gsl_matrix *evec = gsl_matrix_alloc(N,N);
    gsl_vector *noise = gsl_vector_alloc(N);
    
    int iteration=0, it_max=50 ;
    //FILE *ft ; ft = fopen("traj.dat","w");
    while(CONV>CONV_TH && iteration<it_max) {                       // Iterations
        int nh, np ;
        CR_solve(&(D->v),&(D->M),c0);
        
        eigensolve(D->M,eval,evec);             // Noise diagonalization
        init_temp(K);
        for(nh=0;nh<hi0.NH;nh++) {
            for(np=0;np<hi0.NP;np++) {
                
                noise_generation(noise,eval,evec) ;         // Noise generation
                for(int n=0;n<N;n++) D->X.eta[n] = gsl_vector_get(noise,n);
                noise_test(D);                              // Measure noise correlation
                
                traj_solve(D,c0,D->hw[nh].h0) ;
                /*if(np==0) for(int n=0;n<N;n++) { fprintf(ft,"%d\t%lf\t%21.15e\t",iteration,n*p0.dt,D->X.h[n]);
                    for(int m=0;m<N;m+=(int)(N/4)) {
                        if(n>=m) fprintf(ft,"%21.15e\t",D->X.F[tri(n,m)]);
                        else if(n<m) fprintf(ft,"%21.15e\t",1.);
                    }
                    fprintf(ft,"\n");
		    }*/
                measure_kernels(D,K,D->hw[nh].wh0);         // Measure of the new kernels
            }
	    //fprintf(ft,"\n\n");
        }
        //fflush(ft);
        update_kernels(D,K,err,p0);                     // Kernel update
        for(int m=0;m<2;m++) printf("%21.15e\t",err[m]);
        printf("\n\n");
        fflush(stdout);
        CONV = fmax(err[0],err[1]);
        noise_norm(&(D->M),hi0.normNR);                 // Normalize noise correlations
        
        print_output(D->v,D->M,p0,"a") ;                // Output print
	iteration++ ;
    }
    //fclose(ft) ;
    print_output(D->v,D->M,p0,"w") ;                    // Print output functions
    
    gsl_vector_free(eval) ;
    gsl_matrix_free(evec) ;
    gsl_vector_free(noise) ;
     
    return 0 ;
    
}



// FUNCTION DEFINITIONS //


// Initialization of physical parameters, integration constants and initial conditions : from line argument or from code

void init_par(int a, char **argv, phys_t *par, hint_t *hi, ini_t *I0) {
  if(a)
    {
      if(a==1)
	{
	  par->phi = atof(argv[1]) ;          // Rescaled packing fraction
	  par->eps = 1. ;
	  par->b0 = 0. ;
	  par->dt = 0.05 ;
	  par->alpha = 1. ;
	  hi->dh0 = 0.01 ;
	  hi->z = 20. ;
	  hi->threshold = 1e-3 ;
	  hi->NP = 100 ;
	}
      else if(a==9)
        {
	  par->phi = atof(argv[1]) ;          // Rescaled packing fraction
	  par->eps = atof(argv[2]) ;          // Potential energy scale
	  //par->Teff = atof(argv[3]) ;         // Effective temperature
	  //par->taup = atof(argv[4]) ;         // Persistence time
	  par->b0 = atof(argv[3]) ;           // Initial beta
	  par->dt = atof(argv[4]) ;           // Time step
	  par->alpha = atof(argv[5]) ;        // Smooth kernel update
	  //par->N = atoi(argv[5]) ;            // Number of time steps
	  hi->dh0 = atof(argv[5]) ;           // Integration step
	  hi->z = atof(argv[6]) ;             // Rescaled limit of integration
	  hi->threshold = atof(argv[7]) ;     // Error threshold
	  hi->NP = atoi(argv[8]) ;            // Number of realizations
        }
      else
        {
	  printf("Wrong number of arguments!!\n") ;
	  exit(EXIT_FAILURE) ;
        }
    }
  else {
    par->phi = 1.00 ;
    par->eps = 1. ;
    par->b0 = 0. ;
    par->dt = 0.05 ;
    par->alpha = 1. ;
    hi->dh0 = 0.01 ;
    hi->z = 32. ;
    hi->threshold = 1e-3 ;
    hi->NP = 64 ;
  }
  if(par->b0==0){
    I0->k0 = 0. ;
    //I0->MR0 = 0. ;
    //I0->MC0 = 0. ;
    I0->MR0 = par->phi*par->eps*par->eps/2. ;
    I0->MC0 = 2.*I0->MR0 ;
    I0->E0 = 0. ; //par->phi*par->eps/2. ;
    I0->P0 = 0. ; //I0->E0 ;
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
  //hi->hmin = -10. ;
  //hi->hmax = +1. ;
  hi->hmin = -10. ;
  hi->hmax = sqrt(I0->MC0)*par->dt*hi->z + I0->MC0*par->dt*par->dt/2. ;
  hi->NH = (int)((hi->hmax-hi->hmin)/hi->dh0) ;
  hi->norm = par->phi*hi->dh0/(2.*hi->NP) ;
  hi->normNR = 1./((double)(hi->NH*hi->NP)) ;
}


// Initialization of numerical coefficients

void init_coeffs(phys_t par, ini_t I0, coeffs_t *c) {
    c->dt = par.dt ;
    c->dt2 = par.dt*par.dt ;
    c->sM0 = sqrt(I0.MC0) ;
    c->eps = par.eps ;
    c->m0=0 ;
}


// Initialization of dynamical variables

void init_dyn_it(dyn_t *D, phys_t par, hint_t *hi, ini_t I0) {
    
    D->hw = (hweights_t*) calloc(hi->NH,sizeof(hweights_t)) ;
    
    for(int n=0;n<N;n++) {
        D->v.E[n] = 0. ;
        D->v.P[n] = 0. ;
    }
    for(int n=0;n<N1;n++) {
        D->M.R[n]=0. ;
        D->M.C[n]=0. ;
        D->M.chi[n]=0. ;
        D->M.MC[n]=0. ;
        D->M.MCtest[n] = 0. ;
    }

    // Inizialization of the first iteration
    /*
    double A0 = par.phi*par.eps, dtau=par.eps*par.dt ;
    for(int n=0;n<N;n++) {
      for(int m=0;m<=n;m++) {
	D->M.chi[tri(n,m)] = 0.5*A0*(exp(-(n-m)*dtau)-exp(-n*dtau));
	D->M.MC[tri(n,m)] = A0*exp(-(n+m)*dtau) ;
      }
    }
    */
    
    D->M.R[0] = 0.5 ;
    D->M.C[0] = 0.;
    for(int nh=0;nh<hi->NH;nh++) {
        D->hw[nh].h0 = hi->hmin + (nh+0.5)*hi->dh0 ;
        D->hw[nh].wh0 = hi->norm*exp(D->hw[nh].h0-par.eps*bV0(D->hw[nh].h0,par.b0)) ;
    }
    
    print_info(par,*hi,I0) ;                     // Print file with informations
}


void print_info(phys_t par, hint_t hi, ini_t I0) {
 
    FILE *fp ;
    fp = fopen("INFO","w") ;
    
    fprintf(fp,"Physical parameters :\n\n");
    fprintf(fp," phi = %.2e\n eps = %.2e\n\n",par.phi,par.eps) ;
    fprintf(fp," Time step = %.2e\n Steps = %d\n Total time = %.2e\n\n",par.dt,N,(N-1)*par.dt);
    fprintf(fp," Kernel update: alpha = %.2e\n\n",par.alpha);
    fprintf(fp,"Integration over h0 :\n\n") ;
    fprintf(fp," hmin = %lf\n hmax = %lf\n dh = %lf\n NH = %d\n z = %lf\n\n",hi.hmin,hi.hmax,hi.dh0,hi.NH,hi.z) ;
    fprintf(fp," Threshold for the relative error = %.2e\n\n",hi.threshold) ;
    fprintf(fp," Number of realizations for each h0 : NP = %d\n\n",hi.NP) ;
    
    fclose(fp) ;
}

void eigensolve(mat_t M, gsl_vector *eval, gsl_matrix *evec) {
    gsl_matrix *A = gsl_matrix_calloc(N,N) ;
    for(int i=0;i<N;i++) {
           for(int j=0;j<N;j++) gsl_matrix_set(A,i,j,M.MC[tri(i,j)]);
    }
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N) ;
    gsl_eigen_symmv(A,eval,evec,w);
    gsl_eigen_symmv_free(w) ;
}


void noise_generation(gsl_vector *noise,gsl_vector *eval, gsl_matrix *evec){
    
    double threshold=1e-15;
    for(int n=0;n<N;n++) {
        if (fabs(gsl_vector_get(eval,n))<threshold) gsl_vector_set(eval,n,0.) ;
        else if(gsl_vector_get(eval,n)<0) {
            printf("\nWARNING!! Eigenvalue %d = %21.15e\n\n",n,gsl_vector_get(eval,n));
            gsl_vector_set(eval,n,0.) ;
        }
    }
    
    gsl_vector *white = gsl_vector_alloc(N) ;
    for(int n=0;n<N;n++)
        gsl_vector_set(white,n,sqrt(gsl_vector_get(eval,n))*GAUSSGEN) ;
    gsl_blas_dgemv(CblasNoTrans,1.,evec,white,0.,noise);
    //gsl_vector_fprintf(stdout,noise,"%21.15e");
    gsl_vector_free(white);
}


void noise_test(dyn_t *D) {
    for(int n=0;n<N;n++) for(int m=0;m<=n;m++) D->M.MCtest[tri(n,m)] += D->X.eta[n]*D->X.eta[m] ;
}

void noise_norm(mat_t *M, double c) {
    for(int n=0;n<N;n++) for(int m=0;m<=n;m++) M->MCtest[tri(n,m)] *= c;
}



void init_temp(increments_t *K) {
    for(int n=0;n<N;n++) {
        K->E1[n]=0.;
        K->P1[n]=0.;
    }
    for(int n1=0;n1<N1;n1++) {
        K->chi1[n1]=0. ;
        K->MC1[n1]=0. ;
    }
}

void measure_kernels(dyn_t *D, increments_t *K, double c) {
    for(int n=0;n<N;n++) {
        K->E1[n] += c*(D->X.v[n][0]);
        K->P1[n] -= c*(D->X.v[n][1]);
        for(int m=0;m<=n;m++) {
            int nm = tri(n,m) ;
            K->chi1[nm] += c*(D->X.v[n][2]*D->X.F[tri(n,m)]+D->X.v[n][1]) ;
            K->MC1[nm] += c*(D->X.v[n][1]*D->X.v[m][1]) ;
        }
    }
}

void update_kernels(dyn_t *D, increments_t *K, double *err, phys_t par) {
    double /*dk=0., dMR=0.,*/ dMC=0., dchi=0. ;
    for(int n=0;n<N;n++) {
        D->v.E[n] = K->E1[n] ;
        D->v.P[n] = K->P1[n] ;
        for(int m=0;m<=n;m++) {
            int nm = tri(n,m) ;
            dchi += (K->chi1[nm]/D->M.chi[nm]-1.)*(K->chi1[nm]/D->M.chi[nm]-1.) ;
            D->M.chi[nm] = par.alpha*K->chi1[nm] + (1.-par.alpha)*D->M.chi[nm] ;
            dMC += (K->MC1[nm]/D->M.MC[nm]-1.)*(K->MC1[nm]/D->M.MC[nm]-1.) ;
            D->M.MC[nm] = par.alpha*K->MC1[nm] + (1.-par.alpha)*D->M.MC[nm] ;
        }
    }
    dchi /= ((double)N1) ;
    dMC /= ((double)N1) ;

    err[0] = dchi ; err[1] = dMC ;
}


void CR_solve(vec_t *v, mat_t *M, coeffs_t c0) {
 
    int n, l, m ;
    double c1[N]={0.} ;
    for(n=0;n<N;n++) c1[n] = 1./(1.+0.5*c0.dt*M->chi[tri(n,n)]) ;
    
    double r[N1]={0.};
    for(n=0;n<N;n++) {
        int nn = tri(n,n) ;
        r[nn] = -0.5*M->chi[nn] ;
        for(m=n-1;m>=0;m--) {
            double rint = -0.5*M->chi[tri(n,m)]*r[tri(m,m)] ;
            for(l=m+1;l<n;l++) rint -= M->chi[tri(n,l)]*r[tri(l,m)] ;
            r[tri(n,m)] = ( -0.5*M->chi[tri(n,m)] + rint*c0.dt )*c1[n] ;
        }
    }
    for(m=0;m<N;m++) {
        M->R[tri(m,m)] = 0.5 ;
        for(n=m+1;n<N;n++) {
            M->R[tri(n,m)] = M->R[tri(n-1,m)] + 0.5*c0.dt*( r[tri(n-1,m)] + r[tri(n,m)] ) ;
        }
    }
    
    double c[N1]={0.} ;
    c[tri(0,0)] = 0. ;
    for(m=0;m<N;m++) {
        double cint=0.5*(M->MC[tri(0,m)]*M->R[tri(m,m)]+M->MC[tri(0,0)]*M->R[tri(m,0)]) ;
        for(l=1;l<m;l++) cint += M->MC[tri(0,l)]*M->R[tri(m,l)] ;
        c[tri(0,m)] = c0.dt*cint ;
    }
    for(n=1;n<N;n++) {
        for(m=n;m<N;m++) {
            double cint = -0.5*M->chi[tri(n,0)]*c[tri(0,m)] ;
            for(l=m+1;l<n;l++) cint -= M->chi[tri(n,l)]*c[tri(l,m)] ;
            cint += 0.5*(M->MC[tri(n,m)]*M->R[tri(m,m)]+M->MC[tri(n,0)]*M->R[tri(m,0)]);
            for(l=1;l<m;l++) cint += M->MC[tri(n,l)]*M->R[tri(m,l)] ;
            c[tri(n,m)] = c0.dt*cint*c1[n] ;
        }
    }
    for(m=0;m<N;m++) {
        M->C[tri(0,m)] = 0. ;
        for(n=1;n<=m;n++) {
            M->C[tri(n,m)] = M->C[tri(n-1,m)] + 0.5*c0.dt*( c[tri(n-1,m)] + c[tri(n,m)] ) ;
        }
    }
    
    //for(n=1;n<N;n++) v->chiR[n] = v->chiR[n-1] + (M->R[n0]+M->R[n1])*c0.dt/2. ;
}
    
void traj_solve(dyn_t *D, coeffs_t c0, double h0) {
    
    int n, l, m ;
    
    D->X.y[0] = 0. ;
    D->X.h[0] = h0 ;
    for(int m=0;m<3;m++) D->X.v[0][m] = c0.eps*V(h0,m) ; ;
    D->X.F[0] = 1. ;
    
    /*
    // Euler algorithm to compute y(t) and h(t)
    for(n=0;n<N-1;n++) {
        double yint = 0. ;
        for(l=0;l<n;l++) yint -= D->M.chi[tri(n,l)]*(D->X.y[l+1]-D->X.y[l]) ;
        D->X.y[n+1] = D->X.y[n] + c0.dt*( yint - D->X.v[n][1] + D->X.eta[n] ) ;
        D->X.h[n+1] = D->X.h[0] + D->X.y[n+1] + D->M.C[tri(n+1,n+1)] ;
        for(m=0;m<3;m++) D->X.v[n+1][m] = c0.eps*V(D->X.h[n+1],m) ;
    }
    */
    
    double dy[N]={0.} ;
    dy[0] = D->X.eta[0] - D->X.v[0][1] ;
    for(n=0;n<N-1;n++) {
        int n1 = n+1, nn1 = tri(n1,n1) ;
        double ytemp = D->X.y[n] + dy[n]*c0.dt ;
        double hbase = D->X.h[0] + D->M.C[nn1] ;
        double htemp1 = hbase + ytemp ;
        double fint = 0.5*D->M.chi[tri(n1,0)]*dy[0] ;
        for(l=1;l<n1;l++) fint += D->M.chi[tri(n1,l)]*dy[l] ;
        fint = -c0.dt*fint ;
        double c1 = 1./(1.+0.5*c0.dt*D->M.chi[nn1]) ;
        dy[n1] = (fint - c0.eps*V(htemp1,1) + D->X.eta[n1])*c1 ;
        ytemp = D->X.y[n] + 0.5*(dy[n]+dy[n1])*c0.dt ;
        double htemp2 = hbase + ytemp ;
        dy[n1] += c0.eps*c1*(V(htemp1,1)-V(htemp2,1));
        D->X.y[n1] = D->X.y[n] + 0.5*c0.dt*(dy[n]+dy[n1]) ;
        D->X.h[n1] = hbase + D->X.y[n1] ;
        
        for(m=0;m<3;m++) D->X.v[n1][m] = c0.eps*V(D->X.h[n1],m) ;
    }
    
    
    // Evolution of H
    double f[N1]={0.}; double c1[N]={0.};
    for(n=0;n<N;n++) {
        int nn = tri(n,n) ;
        c1[n] = 1./(1.+0.5*c0.dt*(D->M.chi[nn]+D->X.v[n][2])) ;
        f[nn] = -D->X.v[n][2] ;
        for(m=n-1;m>=0;m--) {
            double fint = -0.5*(D->M.chi[tri(n,m)]+D->X.v[n][2])*f[tri(m,m)] ;
            for(l=m+1;l<n;l++) fint -= (D->M.chi[tri(n,l)]+D->X.v[n][2])*f[tri(l,m)] ;
            f[tri(n,m)] = ( c0.dt*fint - D->X.v[n][2])*c1[n] ;
        }
    }
    for(m=0;m<N;m++) {
        D->X.F[tri(m,m)] = 1. ;
        for(n=m+1;n<N;n++) {
            D->X.F[tri(n,m)] = D->X.F[tri(n-1,m)] + 0.5*c0.dt*( f[tri(n-1,m)] + f[tri(n,m)] ) ;
        }
    }
}




// Print of all the outputs files

void print_output(vec_t v, mat_t M, phys_t par, char *mode) {
    
    FILE *fC, *fR, *fD, *fdR, *fk, *fMC, *fMCt, *fchi, *fdMC, *fE, *fP ;
    
    char fileR[FLEN]={'\0'}; filename(fileR,"R",par,mode) ;
    char fileC[FLEN]={'\0'}; filename(fileC,"C",par,mode) ;
    char filedR[FLEN]={'\0'}; filename(filedR,"dR",par,mode) ;
    char fileD[FLEN]={'\0'}; filename(fileD,"D",par,mode) ;
    char filechi[FLEN]={'\0'}; filename(filechi,"chi",par,mode) ;
    char fileMC[FLEN]={'\0'}; filename(fileMC,"MC",par,mode) ;
    char fileMCt[FLEN]={'\0'}; filename(fileMCt,"MCt",par,mode) ;
    char filedMC[FLEN]={'\0'}; filename(filedMC,"dMC",par,mode) ;
    char filek[FLEN]={'\0'}; filename(filek,"k",par,mode) ;
    char fileE[FLEN]={'\0'}; filename(fileE,"E",par,mode) ;
    char fileP[FLEN]={'\0'}; filename(fileP,"P",par,mode) ;
     
    fC = fopen(fileC,mode) ;
    fR = fopen(fileR,mode) ;
    fD = fopen(fileD,mode) ;
    fchi = fopen(filechi,mode) ;
    fMC = fopen(fileMC,mode) ;
    fMCt = fopen(fileMCt,mode) ;
    fdMC = fopen(filedMC,mode) ;
    fk = fopen(filek,mode) ;
    fdR = fopen(filedR,mode) ;
    fE = fopen(fileE,mode) ;
    fP = fopen(fileP,mode) ;
    
    int n,m ;
    for(n=0;n<N;n++) {
        double t = n*par.dt ;
        int nn = tri(n,n) ;
        
        fprintf(fk,"%lf\t%21.15e\n",t,M.chi[nn]);
        fprintf(fE,"%lf\t%21.15e\n",t,v.E[n]);
        fprintf(fP,"%lf\t%21.15e\n",t,v.P[n]);
        fprintf(fdR,"%lf\t%21.15e\n",t,M.R[nn]);
        fprintf(fdMC,"%le\t%21.15e\n",t,M.MC[nn]);
        fprintf(fC,"%lf\t",t);
        fprintf(fR,"%lf\t",t);
        fprintf(fD,"%lf\t",t);
        fprintf(fchi,"%lf\t",t);
        fprintf(fMC,"%lf\t",t);
        fprintf(fMCt,"%lf\t",t);
        
        for(m=0;m<=n;m++) {
            int nm = tri(n,m) ;
            fprintf(fC,"%lf\t",M.C[nm]);
            fprintf(fR,"%lf\t",M.R[nm]);
            fprintf(fD,"%lf\t",M.C[nn]+M.C[tri(m,m)]-2.*M.C[nm]);
            fprintf(fchi,"%lf\t",M.chi[nm]);
            fprintf(fMC,"%lf\t",M.MC[nm]);
            fprintf(fMCt,"%lf\t",M.MCtest[nm]);
        }
        for(m=n+1;m<N;m++) {
            int nm = tri(n,m) ;
            fprintf(fC,"%lf\t",M.C[nm]);
            fprintf(fR,"%lf\t",0.);
            fprintf(fD,"%lf\t",M.C[nn]+M.C[tri(m,m)]-2.*M.C[nm]);
            fprintf(fchi,"%lf\t",M.chi[nn]);
            fprintf(fMC,"%lf\t",M.MC[nm]);
            fprintf(fMCt,"%lf\t",M.MCtest[nm]);
        }
        fprintf(fC,"\n");
        fprintf(fR,"\n");
        fprintf(fD,"\n");
        fprintf(fMC,"\n");
        fprintf(fchi,"\n");
        fprintf(fMCt,"\n");
    }
    fprintf(fk,"\n\n");
    fprintf(fdMC,"\n\n");
    fprintf(fE,"\n\n");
    fprintf(fP,"\n\n");
    fprintf(fC,"\n\n");
    fprintf(fR,"\n\n");
    fprintf(fD,"\n\n");
    fprintf(fchi,"\n\n");
    fprintf(fMC,"\n\n");
    fprintf(fMCt,"\n\n");
    fclose(fk) ; //fclose(fG) ;
    fclose(fdR) ; fclose(fdMC) ;
    fclose(fC) ; fclose(fR) ; fclose(fD) ;
    fclose(fMC) ; fclose(fchi) ;
    fclose(fMCt) ;
    fclose(fE) ; fclose(fP) ;

}


// Standard for output file names

void filename(char *s, char *v, phys_t par,char *mode) {
//sprintf(s,"%s_%.2f_%.0e_%.0e_%.0e_%.0e.dat",v,par.phi,par.eps,par.Teff,par.b0,par.taup);
    if(strcmp(mode,"a")==0) sprintf(s,"%s_%.2f_%.0e_%.0e_it.dat",v,par.phi,par.eps,par.b0);
    else if(strcmp(mode,"w")==0) sprintf(s,"%s_%.2f_%.0e_%.0e.dat",v,par.phi,par.eps,par.b0);
}

