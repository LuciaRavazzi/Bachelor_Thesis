#include <iostream>
#include <cmath>
#include <fstream>
#include "funzioni.h"
#include "nr.h"

			//PER LA PRIMA PARTE DEL PROGRAMMA: RICERCA DELLA FUNZIONE DENSITÀ.


//Variabili e costanti del problema.
  double const a=1.; 			       //Raggio di Bohr [m].	
  double const me=1.;			       //massa elettrone [kg]. 
  double const ht=1.;	    		       //costante di Plank [J*s].
   
  double const A=1; 			       //(2*me)/pow(ht,2);
  double const B=1; 			       //e^2=pow(qe,2)/(4*M_PI*eps);
  double const C=(1./(3*pow(M_PI,2)));

  Dati dati;			
  int N_log = dati.N_log;		       //Punti della mesh logaritmica.
  double Rc = dati.Rc;			       //Raggio di cut-off. 
  int Z = dati.Z;			       //Numero atomico.
			      
  
  Vec_DP X_fissate = mesh(N_log, Rc, dati.Z);  //Punti della mesh totali(logartmica e lineare).

  double rn=(2.36*pow(10,-5)*a)*pow(2*Z,1/3);  //Considero il nucleo esteso: modello a goccia.

  int dim=X_fissate.size();		       //punti della mesh.

  double s=0.;
  double part=0.;
  double sum=0.;
  int k=0; //Indica quante volte ho fatto il ciclo.
  int cost;
  double corrFactor=1.;
  double toll;
  double err=0;
  
  ofstream fileout;
 
//Integrale e integranda.
  funzioneBase *i=new integranda();//definisco l'integranda.
  integrale Hart(X_fissate[0],X_fissate[dim-1],i);//definisco la classe che svolge l'integrale

//Vettori.
  double* nOld=new double[dim];
  double* nNew=new double[dim];
  double* discr=new double[dim];
  double* diff=new double[dim];


//Funzioni.
void Iniz(void);
double potExt(double, int, double rn);
double Particles(double*);
void Print(int);
double Distanza(void);
double Nel(double);
double EKin(void);
double EPot(void);
double PotEl(void);
void PrintEnergy(void);


