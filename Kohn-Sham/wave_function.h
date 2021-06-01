#ifndef wave_function_h
#define wave_function_h

#include "nr.h"
#include "function.h"

class Wave_function:public Function{ //u(r)=R(r)*r;
    
 public:
  Vec_DP Y_rad_D;
    
  Wave_function(const Vec_DP &x, unsigned int ne, unsigned int le, int up_or_down, double Energia=0, unsigned int pieni=0);
  ~Wave_function(){};
		
  double Set_E(double E, ofstream &out3);				//imposta il nuovo autovalore d'energia
  double EigenValue() const;			//restituisce l'autovalore d'energia
		
  double Get_n_pieni()const;		//restituisce il numero di elettroni che occupano la subshell n[l]
  void 	Set_n_pieni(double);    //imposta il numero di elettroni che riempono la subshell n[l]
		
  unsigned int 	Get_l()	const;			//restituisce il momento angolare
  unsigned int	Get_n()	const;			//restituisce il numero quantico principale
  unsigned int	Get_spin() const;		//restituisce lo spin
        
  void Set_l(unsigned int ll);
  void Set_n(unsigned int nn);
        
				
  void Print_EigenValue(const char* file) 	const;	//stampa gli autovalori
  void	Print_EigenFunction(unsigned int n_ciclo, double max)const;	//stampa a ogni ciclo
  void	Print_EigenFunction(double max=0)	const; //stampa finale
		
		
  void operator= (const Wave_function &f);
        
  void Set_Y_rad_D(const Vec_DP &y);
        
  void Print_EigenValue(ofstream &out) const;

  private:
    unsigned int n, l;
    double  n_pieni;				//conservano l'informazione sulla funzione d'onda teorica
    double E_;						//conservano l'autovalore di energia
    int spin; 						//Spin down (up) corrispondono a 0 (1).
};
	
//relazione d'ordine necessaria per riordinare le autofunzioni in funzione della loro energia
bool operator < (const Wave_function &s, const Wave_function &d);
#endif

