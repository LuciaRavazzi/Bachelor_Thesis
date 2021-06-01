
#include "lapack.h"


//Per dsyev
int NB_dsytrd(int N){
    
    int ispec= 1;
    char name[7] = "dsytrd";
    char opts[2] = "U";
    int n2=-1, n3=-1, n4=-1; //Parametri fissati
    
    //Blocksize 
    int NB_dsytrd = ilaenv_(&ispec,name,opts,&N,&n2,&n3,&n4);
    
    return NB_dsytrd;

}


//Per dsysv
int NB_dsytrf(int N){
    
    int ispec= 1;
    char name[7] = "dsytrf";
    char opts[2] = "U";
    int n2=-1, n3=-1, n4=-1; //Parametri fissati
    
    //Blocksize
    int NB_dsytrf = ilaenv_(&ispec,name,opts,&N,&n2,&n3,&n4);
    
    return NB_dsytrf;
    
}




//Output
//Autovalori in ordine crescente nel vettore w
//Autovettori per righe nella matrice a (SOVRASCRITTA). Sono nell'ordine degli autovalori

void diagonalize(int &N, double *a, double *w, int &lwork, double *work, int &info) {

    char jobz= 'V';
    char uplo= 'U';
    int lda= N;

    dsyev_(&jobz, &uplo, &N, a, &lda, w, work, &lwork, &info);

}



//Output
//In b è salvata la matrice inversa di a (quindi b è SOVRASCRITTA). Attenzione: anche a viene sovrascritta

void inverse(int &N, double *a, int *ipiv, double *b, int &lwork, double *work, int &info){
   
    char uplo= 'U';
    int nrhs= N;
    int lda= N;
    int ldb= N;
    
    //Sovrascrivo b con la matrice identità
    for(int i=0; i<N; ++i)
        for(int j=0; j<N; ++j)
            if (j==i)
                b[i*N+j]=1;
            else
                b[i*N+j]=0;
    
    dsysv_(&uplo,&N,&nrhs,a,&lda,ipiv,b,&ldb,work,&lwork,&info);

}


