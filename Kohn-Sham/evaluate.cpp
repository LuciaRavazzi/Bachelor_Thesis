#include "V_Diretto.h"
#include "evaluate.h"

extern int N;
extern Dati dati;
extern double *rad_D;
extern double *rad_inv_D;
extern double *G_inv;

extern int lwork;
extern double *work;
extern int info;


extern double *matrice;
extern double *autovalori;

extern double *T;
extern double *T_prova;

extern unsigned int lmax;
extern unsigned int nmax;

extern unsigned int num_cicli;



double Evaluate_wave_function(Density &d, Potenziale &potenziale, ofstream &out2, ofstream &out3){

 unsigned int cont=0;  //contatore delle funzioni d'onda corrette per l fissato e spin fissato
 double err=0;	//tiene conto della variazione quadratica dell'autovalore	
 ofstream out_matrice;
           
 for(unsigned int l=0; l<=lmax; ++l)	{		//ciclo sul numero quantico l   
  for(int spin=0; spin<=1; ++spin){			//loop sullo spin          
    //Risetto gli elementi di matrice a zero
    for(int i=0; i<N; ++i)
     for(int j=0; j<N; ++j)
       matrice[ind(i,j)]=0;     
       //Potrei dimezzare le operazioni perchè matrice è simmetrica
        for(int i=0; i<N; ++i){
         matrice[ind(i,i)]= l*(l+1)/(2.*pow(d.X[i],2))+potenziale.Y[i];
            for(int j=0; j<N; ++j)
                matrice[ind(i,j)]+= T[ind(i,j)];
         }
         
 
   // for(int j=0; j<N; ++j){
   //  cout << potenziale.Y[i] << " ";   
   // }
   
//Al primo ciclo di HF le autofunzioni sono identicamente nulle. Quindi rimane l'equazione agli autovalori per l'atomo a un elettrone con nucleo esteso di carica Z
        diagonalize(N, matrice, autovalori, lwork, work, info);
        cout << "INFO " << info << endl;   
                                            
//Vec_DP di appoggio per copiare i valori assunti dalla singola autofunzione nei punti della mesh
   Vec_DP valori_Y_rad_D(N);
   Vec_DP valori_Y(N);
        
   double a; //per normalizzazione
   double incremento_errore; //incremento errore su autovalori
                           
   for(int n=l+1; n<=nmax; ++n){            
    for(int j=0; j<N; ++j)
          	valori_Y_rad_D[j]= matrice[ind(cont,j)];
     a= copysign(1.,valori_Y_rad_D[0])/sqrt(norm2(valori_Y_rad_D));
            
    for(int j=0; j<N; ++j)
     valori_Y_rad_D[j]*=a;
            
    for(int j=0; j<N; ++j)
     valori_Y[j]= rad_inv_D[j]*valori_Y_rad_D[j];
                
    //Scelgo la funzione d'onda corrispondente del ciclo precedente e la aggiorno
    for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k)
     if(d.psi[k].Get_l()==l && d.psi[k].Get_spin()==spin && d.psi[k].Get_n()==n){
      //Valori della funzione d'onda nella mesh
      d.psi[k].Set_Y_rad_D(valori_Y_rad_D);
      d.psi[k].Set_Y(valori_Y); //R(r)*r
      //Interpolo
      d.psi[k].Set();
      //cout << n << " " << l << " " << spin << " autovalori[cont] " << autovalori[cont] << endl;     
      incremento_errore= d.psi[k].Set_E(autovalori[cont], out3);
      out2 << l << " " << spin-0.5 << " " << n << " " << incremento_errore <<"      ";
    
      if(d.psi[k].Get_n_pieni()!=0)
      err+= incremento_errore;
     }
    cont++;
   }//loop su n
   cont=0;
  }//loop sullo spin
 } //loop su l   
    
    
 cout <<  "lmax " << lmax << endl;
 cout << "nmax " << nmax << endl;              
 d.Evaluate_Y(); //aggiorno density e riordino le funzioni d'onda
 return sqrt(err/d.Get_n_autofunzioni());                  
}
               

double norm2(Vec_DP v){
	double sum=0;
	for(int i=0; i<v.size(); i++)
    	sum+=v[i]*v[i];
	return sum;
}
             
