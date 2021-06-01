
#ifndef lapack_h
#define lapack_h

#include "nr.h"


extern "C" {
    int dsyev_(char *JOBZ,
               char *UPLO,
               int *N,
               double *A,
               int *LDA,
               double *W,
               double *WORK,
               int *LWORK,
               int *INFO);
    
    
    
    int dsysv_(char *UPLO,
               int *N,
               int *NRHS,
               double *A,
               int *LDA,
               int *IPIV,
               double *B,
               int *LDB,
               double *WORK,
               int *LWORK,
               int *INFO);
    
        
    
    int ilaenv_(int *ISPEC,
                char *NAME,
                char *OPTS,
                int *N1,
                int *N2,
                int *N3,
                int *N4);
    
}


int NB_dsytrd(int N);

int NB_dsytrf(int N);



void diagonalize(int &N, double *a, double *w, int &lwork, double *work, int &info);



void inverse(int &N, double *a, int *ipiv, double *b, int &lwork, double *work, int &info);




#endif



