/* header for pseudorandom numbers generator */

/*********************************************************/
/* UNIGEN is the uniform pseudorandom numbers generator  */
/*         can opt changing the GENSET definition        */
/* GAUSSGEN is the Gaussian pseudorandom generator       */
/* INIRAND is the seed initializer and must be called    */
/*         at the beginning                              */
/*********************************************************/

/* possible values for GASGEN: GRAND */
#define GAUSSGEN GRAND    /* gaussian random generator */

/* possible values for GENSET: MYRANF, DRAND, MERSENNE */
#define GENSET MERSENNE   /* uniform random generator  */

//#define SEED 11061973

/******* DEFINE THE FOLLOWING FOR MERSENNE WITH TIME INITIALIZATION */
#define TIMESEED 

/**********************************************************************/

#define GRAND grand(0.0,1.0)

#define MYRANF 1
#define DRAND 2
#define MERSENNE 3


/* seed and generator initialization */
#if GENSET==MYRANF 
#define INIRAND myseed
#define UNIGEN myranf()
#elif GENSET==DRAND
#define INIRAND srand48
#define UNIGEN drand48()
#elif GENSET==MERSENNE
#define INIRAND sgenrand
#define UNIGEN genrand()
#endif


/************************** MERSENNE PARAMETERS ********************/
/* Period parameters */  
#define N_MERS 624
#define M_MERS 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

/********************************************************************/


double myranf(void);    /* dumb generator */
void myseed(unsigned long seed);    /* dumb initializer */


void sgenrand(unsigned long seed);  /* mersenne initializer */ 
double genrand(void);    /* mersenne twister */

double grand(double m, double s);    /* gaussian generator */
