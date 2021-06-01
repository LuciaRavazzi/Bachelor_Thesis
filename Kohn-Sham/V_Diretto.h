#ifndef V_Diretto_h
#define V_Diretto_h

#include "nr.h"
#include "function.h"
#include "density.h"

using namespace std;

//Tiene conto del potenziale di Hartree.

class V_Diretto: public Function {
 public:
  //Costruttore di potenziale nullo.
  V_Diretto(const Vec_DP &x);   					    
  //Aggiorno dopo aver ricevuto la densità.
  V_Diretto(const Density &);  
  //Altro costrutture.
  V_Diretto(const Vec_DP &x, const Vec_DP &y, double yp1=0, double ypn=0);  	     
  //Valuta la funzione potenziale nei punti della mesh.      
  void Evaluate(const Density &);	
  //Uguaglia due funzioni.						
  void operator= (const V_Diretto &f);
};

V_Diretto operator*(const double, V_Diretto &);
V_Diretto operator+(V_Diretto , V_Diretto );
V_Diretto operator-(V_Diretto , V_Diretto );

//restituisce la differenza fra 2 potenziali
double Get_err(V_Diretto &s, V_Diretto &d); 


#endif

