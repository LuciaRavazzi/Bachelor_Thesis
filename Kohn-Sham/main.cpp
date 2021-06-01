#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "nr.h"
#include "lapack.h"
#include "funzioni.h"
#include <algorithm> 
#include <vector>
#include <iomanip>
#include <string>
#include "density.h" 
#include "potential.h"
#include "e_var.h"
#include "evaluate.h"
#include "V_Diretto.h"

//0.000227
#define R0  0.0000236 

using namespace std;


//Memorizzo i parametri passati da file
Dati dati;			

int N_log = dati.N_log;

//Raggio di cutoff
double Rc = dati.Rc;

Vec_DP X_fissate = mesh(N_log, Rc, dati.Z);


//Variabili extern in alcune librerie
int N = X_fissate.size(); //ordine delle matrici usate, pari al numero di punti della mesh
double *G_inv;
double *D;
double *rad_D;
double *rad_inv_D;
double *fact;

double *matrice;
double *autovalori;

double *C;
double *T;
double *T_prova;

unsigned int lmax;
unsigned int nmax;

unsigned int num_cicli;

//Potenziale del nucleo
Potenziale V_0(X_fissate);


// Dichiarazione delle variabili per il funzionamento di LAPACK. Usate come variabili extern in alcune librerie

//Valori di lwork per efficienza ottimale. È lecito scegliere valori minori per risparmiare memoria, ma l'efficienza non è più ottimale. In ogni caso si deve avere lwork_dsyev>=(3N-1) e lwork_dsysv>=1
int lwork_dsyev = (NB_dsytrd(N)+2)*N;
int lwork_dsysv = (NB_dsytrf(N))*N;

int lwork = max(lwork_dsyev,lwork_dsysv);

double *work = new double[lwork];

//Per dsysv
int *ipiv = new int[N];

int info;



int main() {

 cout << "Aperto " << endl;
 double tempo_tot=0;	//tiene il tempo totale del programma in sec.
 double tempo;	//per intervalli parziali
 clock_t start,end;
 start=clock();

 cout << "Z= " << dati.Z << " N= " << dati.N << endl;	


//Calcolo i fattoriali per i coefficienti di CG
 fact= new double[31];
 set_factorial(fact);

        
//Scrivo la matrice dei pesi D, D^0.5, D^-0.5 (matrici diagonali).
//Metodo dei trapezi    
 D= new double[N];
 rad_D= new double[N];
 rad_inv_D= new double[N];
    
 D[0]= X_fissate[1]/2.;
 D[N-1]= (Rc-X_fissate[N-2])/2.;
    
 for(int i=1; i<N-1; ++i)
  D[i]= (X_fissate[i+1]-X_fissate[i-1])/2.;
    
 for(int i=0; i<N; ++i){
  rad_D[i] = sqrt(D[i]);    
  rad_inv_D[i] = 1./rad_D[i];
 }
        
//Matrice C tridiagonale NON simmetrica che approssima derivata seconda
        
 C= new double[N*N];
    
 C[ind(0,0)]= -2./((X_fissate[1]-X_fissate[0])*X_fissate[0]);
 C[ind(N-1,N-1)]= -2./((Rc-X_fissate[N-1])*(X_fissate[N-1]-X_fissate[N-2]));
    
 C[ind(0,1)]= 2./((X_fissate[1]-X_fissate[0])*X_fissate[1]);
 C[ind(1,0)]= 2./((X_fissate[1]-X_fissate[0])*(X_fissate[2]-X_fissate[0]));
                    
 C[ind(N-2,N-1)]= 2./((X_fissate[N-1]-X_fissate[N-2])*(X_fissate[N-1]-X_fissate[N-3]));
    
 C[ind(N-1,N-2)]= 2./((X_fissate[N-1]-X_fissate[N-2])*(Rc-X_fissate[N-2]));
             
             
 for(int i=1; i<=N-2; ++i){
  C[ind(i,i)]= -2./((X_fissate[i+1]-X_fissate[i])*(X_fissate[i]-X_fissate[i-1]));
 }
    
 for(int i=1; i<=N-3; ++i){
  C[ind(i,i+1)]= 2./((X_fissate[i+1]-X_fissate[i])*(X_fissate[i+1]-X_fissate[i-1]));
  C[ind(i+1,i)]= 2./((X_fissate[i+1]-X_fissate[i])*(X_fissate[i+2]-X_fissate[i]));        
 }
 
        
//Matrice cinetica con differenze finite e metodo dei trapezi
 T= new double[N*N];
    
 T[ind(0,0)]= 1./((X_fissate[1]-X_fissate[0])*X_fissate[0]);
 T[ind(N-1,N-1)]= 1./((Rc-X_fissate[N-1])*(X_fissate[N-1]-X_fissate[N-2]));
    
 T[ind(0,1)]=T[ind(1,0)]= -1./((X_fissate[1]-X_fissate[0])*
    				sqrt(X_fissate[1]*(X_fissate[2]-X_fissate[0])));
                    
 T[ind(N-2,N-1)]=T[ind(N-1,N-2)]= -1./((X_fissate[N-1]-X_fissate[N-2])*
    				sqrt((X_fissate[N-1]-X_fissate[N-3])*(Rc-X_fissate[N-2])));
                          
 for(int i=1; i<=N-2; ++i){
  T[ind(i,i)]= 1./((X_fissate[i+1]-X_fissate[i])*(X_fissate[i]-X_fissate[i-1]));
 }
    
 for(int i=1; i<=N-3; ++i){
  T[ind(i,i+1)]=T[ind(i+1,i)]= -1./((X_fissate[i+1]-X_fissate[i])*
 				sqrt((X_fissate[i+1]-X_fissate[i-1])*(X_fissate[i+2]-X_fissate[i])));
 }
    
       

//Matrice e vettore per diagonalize, definiti extern in density.cpp
 matrice = new double[N*N];
 autovalori = new double[N];
//Inizializzo a 0 gli elementi di matrice
 for(int i=0; i<N; ++i)
   for(int j=0; j<N; ++j)
     matrice[ind(i,j)]=0;
           
 create_folder("risultati");
 create_folder("risultati/autofunzioni");
      
    
//Densità: all'inizio è nulla ovunque.
 Density density(X_fissate, dati.Z, dati.N);
 //Density density_new(X_fissate, dati.Z, dati.N);
 //density.Evaluate_Y_int(); //Se volessi partire da una funzione densità esponenziale.
 

//Calcolo il potenziale del nucleo nei punti della mesh
 Vec_DP nucleo(N);
//Stima del raggio del nucleo
 double Rn= R0*pow(2*density.Get_Z(),(double)1./3.);
//Calcolo dell'energia potenziale di interazione con il nucleo
 for(int i=0; i<X_fissate.size(); ++i){
  if(X_fissate[i]<=Rn)
   nucleo[i]=-(double)1.5*density.Get_Z()/Rn+(double)0.5*density.Get_Z()*pow(X_fissate[i],2)/pow(Rn,3);
  else
   nucleo[i]=-(double)density.Get_Z()/X_fissate[i];
 }

//Viene valutato direttamente dal costruttore e aggiornato a ogni iterazione di HF (con un mixing)
 Potenziale potenziale(density, nucleo); //Per adesso ho solo il potenziale del nucleo.

 V_Diretto V_diretto(density);  //Nullo.

 lmax= density.Get_lmax();
 nmax= density.Get_nmax(); 

//Salvo il potenziale iniziale, cioè il potenziale del nucleo perchè la densità è nulla.
 V_0=potenziale;
 V_0.Print("risultati/potenziale0.dat");

//Codice che implementa il ciclo di HF	
 int controllo=1;	
 double controllo_err;		
 num_cicli=0;
 double err_pot;
 double E_var=0,E_var_old,err_E_var,E_var_appo;
 unsigned int contatore= 1; //conta il numero di volte che ricomincio a ciclare con altre 50 iterazioni


 end=clock();
 tempo_tot+=((double)(end-start))/CLOCKS_PER_SEC;
	
 ofstream out;
 out.open("elenco.dat");
    
 ofstream out2;
 out2.open("elenco2.dat");
    
 ofstream out3;
 out3.open("elenco3.dat");

//Nello step 0 valuto autovalori e autofunzioni di prova con un potenziale quasi coulombiano (nucleo esteso)
 cout << "Inizio il " << num_cicli << " ciclo di HF" << endl;
    
//Stampa numero ciclo
 out << num_cicli<< "     " << setprecision(7);
 out2 << num_cicli<< "     " << setprecision(7);
 out3 << num_cicli<< "     " << setprecision(7);

//Calcolo autofunzioni e autovalori di prova (nucleo esteso)
 controllo_err= Evaluate_wave_function(density, potenziale, out2, out3);
 cout << "Errore al primo ciclo " << controllo_err << endl;

//Calcolo energia variazionale (somma degli autovalori)
  //E_var_appo= Evaluate_E_var1(density,V_diretto,scambio,out);    
  E_var= Evaluate_E_var2(density,V_diretto,out,0,dati.Z,dati.N); //Energia variazionale.
  cout << "Energia variazionale al primo ciclo " << E_var << endl;

//Stampa autofunzioni e autovalori di prova
  //density.Print_EigenFunction(num_cicli);
  density.Print_EigenValue(out);
  out << endl << endl;
  out2 << endl << endl;
  out3 << endl << endl; 

//Potenziali per il mixing.
 Potenziale potenziale_new(density,nucleo);

 V_Diretto V_diretto_new(density);

  
 do {  
   num_cicli++;
        
   start=clock();
   cout << "Inizio il " << num_cicli << " ciclo di HF" << endl;
                
   //Stampa numero ciclo
   out << num_cicli<< "     ";
   out2 << num_cicli<< "     ";
   out3 << num_cicli<< "     ";
         
   //Aggiorno il potenziale efficace 
   potenziale_new.Evaluate(density,nucleo);   
   V_diretto_new.Evaluate(density);	 	
  

   //Mixing
   potenziale= (1-dati.alpha)*potenziale+dati.alpha*potenziale_new;
   V_diretto=(1-dati.alpha)*V_diretto+dati.alpha*V_diretto_new;
   

   //Density viene aggiornata con nuovi autovalori e autofunzioni
   controllo_err= Evaluate_wave_function(density, potenziale, out2, out3);

   if(info!=0)
     cout<< "Errore"<< "    " <<"INFO= " << info << endl;
     out2 << endl << endl;
     out3 << endl << endl;
 
    //Calcolo di E_var con potenziale mixati
    E_var_old=E_var;	//Salvo il valore precedente di energia variazionale

    E_var= Evaluate_E_var2(density,V_diretto,out,0,dati.Z,dati.N);   

    err_E_var=fabs(E_var-E_var_old);

        
    //Stampa autofunzioni e autovalori
    //density.Print_EigenFunction(num_cicli);
    density.Print_EigenValue(out);
    out << endl << endl;
        
            
    //Da sistemare con anche termine di Fock??
    err_pot=Get_err(potenziale, potenziale_new); //Ne estraggo l'errore
        
        
    end=clock();
    tempo=((double)(end-start))/CLOCKS_PER_SEC;  //Valuto il tempo dello step
		
    cout << setprecision(10) << "E_totale= "<<E_var<< endl <<"err_E_totale= "<<err_E_var <<endl << "tempo dello step= "<<tempo<<endl;
    cout<<"err_pot= "<<err_pot<< endl << "err_autoval= "<<controllo_err<<endl<<endl;
    tempo_tot+=tempo;
    if (num_cicli==dati.controllo_cicli*contatore) {
     cerr<<" continuare a ciclare? (1/0)"<<endl;
     cin>>controllo;  
     if(controllo==1)
      contatore++;
     }

 }while( (err_E_var>dati.err || controllo_err>dati.err) && controllo==1);
	
 E_var= Evaluate_E_var2(density,V_diretto,out,1,dati.Z,dati.N);  //Valuto l'energia alla fine del ciclo e scrivo su file.
 start=clock();
	
 //Stampa a completamento del programma
 density.Print("risultati/densita_f.dat");
 potenziale.Print("risultati/potenziale_f.dat");
 density.Print_EigenValue("risultati/autovalori_f.dat");
 density.Print_EigenFunction();
    
    	
 //Stampa dei risultati globali
 ofstream out_finale;
 out_finale.open("risultati/energia.dat");
 out_finale<< setprecision(12) <<"E totale= "<<E_var<<endl;
 out_finale<< setprecision(3) <<"err E tot= " <<err_E_var<< endl << "err pot= "<<err_pot<< endl <<"err_autoval= "<<controllo_err<<endl;
	
 end=clock();
 tempo_tot+=((double)(end-start))/CLOCKS_PER_SEC;
	
 out_finale<<"L'algoritmo ha compiuto " << num_cicli << " iterazioni di HF in " <<tempo_tot <<" s"<<endl;
 out_finale.close();

 ofstream fileout;
 fileout.open("risultati/Energia_variazionali_finali.dat",ios::app);
 fileout << dati.Z << " " << dati.N << " " << E_var << " " << err_E_var << endl;
 fileout.close();

    return 0;
}






