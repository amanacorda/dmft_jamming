/* header for random numbers generation */
#include<math.h>
#include<stdlib.h>
#include "random.h"
#include<time.h>

static unsigned long iseed;     /* dumb seed */
static unsigned long mt[N_MERS]; /* the array for the state vector  */
static int mti=N_MERS+1; /* mti==N_MERS+1 means mt[N_MERS] is not initialized */
/******************************** M Y R A N F **********************/

void
myseed(unsigned long seed)
{
  iseed=seed;
}


double myranf(void)
{
  iseed*=65549;
/*    if (iseed<0) iseed+=2147483647; */
  return .465661288e-9*iseed;
}
  

/********************** M E R S E N N E   T W I S T E R  ***************/
/* Initializing the array with a seed */
void sgenrand(unsigned long seed) 
{
    int i;

    for (i=0;i<N_MERS;i++) {
         mt[i] = seed & 0xffff0000;
         seed = 69069 * seed + 1;
         mt[i] |= (seed & 0xffff0000) >> 16;
         seed = 69069 * seed + 1;
    }
    mti = N_MERS;
}

 /* generating reals */
double genrand(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N_MERS) { /* generate N_MERS words at one time */
        int kk;

        if (mti == N_MERS+1)   /* if sgenrand() has not been called, */
#ifdef TIMESEED
	  sgenrand((unsigned long)time(0));
#else
	sgenrand(4357); /* a default initial seed is used   */
#endif

        for (kk=0;kk<N_MERS-M_MERS;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M_MERS] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N_MERS-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M_MERS-N_MERS)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N_MERS-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N_MERS-1] = mt[M_MERS-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y * 2.3283064365386963e-10 ); /* reals: [0,1)-interval */
    /* return y; */ /* for integer generation */
}

/* This main() outputs first 1000 generated numbers.  */
/* main() */
/* {  */
/*     int i; */

/*     sgenrand(4357); */
/*     for (i=0; i<1000; i++) { */
/*         printf("%10.8f ", genrand()); */
/*         if (i%5==4) printf("\n"); */
/*     } */
/* } */



/************************** G R A N D *****************************/

double 
grand(double m, double s) /* normal random variate generator */ 
{ /* mean m, standard deviation s */
  double x1, x2, w, y1; 
  static double y2; 
  static int use_last = 0; 
 
  if (use_last) /* use value from previous call */
    {
      y1 = y2; 
      use_last = 0; 
    }
  else
    {
      do {
	x1 = 2.0 * UNIGEN - 1.0; 
	x2 = 2.0 * UNIGEN - 1.0; 
	w = x1 * x1 + x2 * x2; 
      } while ( w >= 1.0 ); 

      w = sqrt( (-2.0 * log( w ) ) / w ); 
      y1 = x1 * w; 
      y2 = x2 * w; 
      use_last = 1; 
    }

  return( m + y1 * s );  
}



