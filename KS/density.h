#ifndef density_h
#define density_h

/*
La classe Density al suo interno contiene un vettore con gli orbitali pieni.
Ne gestisce in automatico aggiornamento e somma in modo tale da avere un'unica  
funzione densità (interpolata)
*/

#include "nr.h"
#include "lapack.h"
#include "funzioni.h"
#include "function.h"
#include "wave_function.h"
#include "my_nr.h"


#include <vector>
using namespace std;

class Density: public Function
	{
	private:
 
        unsigned int n_autofunzioni; 						//Quanti orbitali sto studiando?       		
        void Put_e();								//Ordina gli orbitali e li riempe
	void Def_wave_function();						//Mi definisce gli orbitali
	unsigned int 	Get_periodo() const;					//Mi dice fino a che periodo devo creare gli orbitali
	
        
        unsigned int 	Z_;										//contiene il numero totale di elettroni
	unsigned int 	N_;										//contiene il numero totale di elettroni
		
        
        unsigned int nmax; //massimo valore del numero quantico principale
        unsigned int lmax; //massimo valore del numero quantico l
        
        
        
 public:
	Density(const Vec_DP &x, unsigned int Z, unsigned int N);
        Density(const Vec_DP &x, const Vec_DP &y, double yp1=0, double ypn=0); 
	~Density(){};
		
     
        vector<Wave_function> psi; //Vettore che raccoglie le funzioni d'onda di ciascuna subshell

	void 	Evaluate_Y(); //Riassume i risultati valutando la distribuzione totale
        void Evaluate_Y_int();	
		
        void 	Print_EigenValue(const char* file) const; //Stampa gli autovalori in modo carino
	void 	Print_EigenFunction(unsigned int n_ciclo, double max=0) const; //Stampa le autofunzioni a ogni ciclo 
	void	Print_EigenFunction(double max=0) const; //Stampa autofunzioni finali
        
	unsigned int Get_Z() const; //Restituisce il numero di protoni (controllo che il nucleo si integri bene)
	unsigned int Get_N() const; //Restituisce il numero di elettroni schermanti (controllo che la schermatura si integri bene)
        unsigned int Get_n_autofunzioni() const;
		
	double Integrale(const unsigned int n) const; //Restituisce la primitiva con costante iniziale 0. 
		
        void operator= (const Density &f);

        unsigned int Get_nmax();
        unsigned int Get_lmax();
        
        
        void Print_EigenValue(ofstream &out) const;
		
	};
    
    
Density operator*(const double val, Density &f);
Density operator+(Density s, Density d);
Density operator-(Density s, Density d);
   
#endif

