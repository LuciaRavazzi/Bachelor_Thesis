
#include "integrale.h"
#include <iostream>
#include <cmath>


using namespace std;

integrale::integrale(double a,double b,funzioneBase *f){
	if(a<b){
		_a=a;
		_b=b;
	} else {
		_a=b;
		_b=a;
	}	
		
	_f=f;
}


double integrale::trapezi(Vec_DP& passo, double* n, int dim, double r){
//double* passo definisce la mesh.
//double* n è il vettore con le componenti della densità.
//dim è il numero di punti della mesh.
//r è la distanza rispetto al quale voglio calcolare la funzione densità (è fissata quando faccio l'integrale).

	double sum=0.;

	for(int i=0; i<dim-1; i++){ //dim-1 per _f->eval(passo[i+1],n).
		if(r>=passo[i]){
	         sum+=(1./r)*(_f->eval(passo[i+1],n[i+1],2)+_f->eval(passo[i],n[i],2))*(passo[i+1]-passo[i]);
		} else { 
                 sum+=(_f->eval(passo[i+1],n[i+1],1)+_f->eval(passo[i],n[i],1))*(passo[i+1]-passo[i]);
		}
        } 
		
       return (4*M_PI)*sum/2.;
}
	





