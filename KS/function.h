#ifndef function_h
#define function_h

#include "nr.h"


//Potrei ripristinare i data membri protected o forse anche private


//Posso definire static Vec_DP X e risparmiare memoria
class Function
	{
	public://per arrvarci anche nelle classi derivate
		Vec_DP X;		Vec_DP Y;		Vec_DP Y2;		Vec_DP Yp;
	
        
		Function(const Vec_DP &x);
		Function(const Vec_DP &x, const Vec_DP &y, double yp1=0, double ypn=0);
		Function(const Function &f);
		~Function(){};
			
		//Restituisce il vettore delle ascisse
		Vec_DP	Get_X() const;
		//Vec_DP	Get_Y() const;
		
		//Svolge lo spline a partire dai dati inseriti
		void 	Set_Y(const Vec_DP &y, double yp1=0, double ypn=0);
		void	Set(double yp1=0, double ypn=0);
		
		//Restituisce il numero di punti su cui si sviluppa l'interpolazione
		unsigned int size() const;
		
		//Funzioni per la stampa
        void	Print(const char*) const;
		void	Print(const char*, double, double =0.001) const;
			
		//Funzione per l'estrazione delle valutazioni
		virtual double 	Value(const double x) const;		//Restituisce la valutazione della funzione in x
		double 	D1(const double x) const;					//Permette di valutare la derivata prima (non utilizzato)
		double 	D2(const double x) const;					//Permette di valutare la derivata seconda (non utilizzato)
		double	Primitive(const unsigned int n) const;		//Permette di valutare la primitiva con costante iniziale 0
		double 	Normalize();								//Normalizza le funzioni
		
		double operator() (const double x) const;			//Forma compatta di value
		double& operator[] (const unsigned int i);			//Referenza all'i-esimo punto valutato
		void operator= (Function &f);
	};

Function operator*(const double, Function &);
Function operator+(Function &, Function &);
Function operator-(Function &, Function &);

#endif

