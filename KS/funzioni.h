#ifndef funzioni_h
#define funzioni_h

#include "nr.h"

extern int N;



//Struttura di comodo per gestire l'ingresso di dati
class Dati
{
  public:
  	Dati();
    unsigned int N_log;	
    double Rc;
    unsigned int Z;
    unsigned int N;	
    
    double err; 
    double accuratezza;
    unsigned int controllo_cicli;
    double alpha;
    
    //void load(); 	//Costruttore a partire dal file
};



Vec_DP mesh(int N_log, double Max, unsigned int Z);

Vec_DP mesh_lin(int N, double Max, unsigned int Z);

Vec_DP mesh_log(int N, double Max, unsigned int Z);



//Restituisce la posizione dell'elemento di matrice (i,j) all'interno dell'array. Nel mio codice le matrici sono trattate come semplici array per avere compatibilità con LAPACK. N.B: i e j da 0 a (N-1)
inline unsigned int ind(unsigned int i, unsigned int j) {
    return i*N+j;
}

//Restituisce la posizione dell'elemento di matrice (i,j) all'interno di un vector di oggetti di tipo Exchange.
inline unsigned int ind2(unsigned int i, unsigned int j) {
    return i*2+j;
}



//Calcolo dei fattoriali
double factorial(unsigned int k);


//Conviene calcolare una volta per tutti i fattoriali che possono servire a valutare i coefficienti di CG
void set_factorial(double *fact);

//Per creare una cartella
void create_folder(char *nomefile);


#endif
