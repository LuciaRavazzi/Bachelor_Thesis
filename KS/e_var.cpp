
#include "e_var.h"

extern double *T;
extern double *D;
extern double *rad_D;
extern double *rad_inv_D;
extern double *G_inv;
extern Dati dati;
extern Potenziale V_0;
extern double *T_prova;
extern double *C;


double Evaluate_E_var1(Density &d, Potenziale &V,  ofstream &out){
    
    //Somma degli autovalori
    double sum_autovalori=0;
	for (unsigned int i=0; i<d.Get_n_autofunzioni(); ++i) 
    	sum_autovalori+= d.psi[i].Get_n_pieni()*d.psi[i].EigenValue();
    //cout << "Somma autovalori = "<< sum_autovalori << endl;
    //out << sum_autovalori << "     ";
    double E_var= sum_autovalori;
    
    //Termine diretto
    double sum_diretto=0;
    for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k){
      
      //Non necessario ma velocizza
      if (d.psi[k].Get_n_pieni()>0){
      
        //Somme su i
        for (unsigned int i=0; i<d.size(); ++i)
        	sum_diretto+= -0.5*d.psi[k].Get_n_pieni()*V[i]*d.psi[k].Y_rad_D[i]*d.psi[k].Y_rad_D[i];
       }
     }
	//cout << "Termine diretto = "<< sum_diretto << endl;
    //out << sum_diretto << "     ";
    E_var+=sum_diretto;
    
/*    //Termine derivante dal potenziale di scambio
    double sum_scambio=0;
    unsigned int l;
    int spin;

	for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k){
      //Non necessario ma velocizza
      if (d.psi[k].Get_n_pieni()>0){
      
    	//Valore di l
        l= d.psi[k].Get_l();
        spin= d.psi[k].Get_spin();
    	        
        //Somme su i e j
        for (unsigned int i=0; i<d.size(); ++i)
        	 for (unsigned int j=0; j<d.size(); ++j)
        sum_scambio+=-0.5*d.psi[k].Get_n_pieni()* 
                    s[ind2(l,spin)].K[ind(i,j)]*d.psi[k].Y_rad_D[i]*rad_D[i]*
                    		d.psi[k].Y_rad_D[j]*rad_D[j]; 
      }
     }
    //cout << "Termine Fock = "<< sum_scambio << endl; 
    //out << sum_scambio << "     ";
    E_var+=sum_scambio;
*/
    out << E_var << "                    ";    
    
	return E_var;
    
}




double Evaluate_E_var2(Density &d,V_Diretto &V_diretto, ofstream &out, int condition, int Z, int N){

 double E_var=0;
 double sum_autovalori=0;  

//Valuto l'occupazione delle autofunzioni.
 
 for(int i=0; i<d.Get_n_autofunzioni(); i++){
  double l=d.psi[i].Get_l();
  double n=d.psi[i].Get_n();
  char a,b;
  if(l == 0)
    a = 's';
  if(l == 1)
    a = 'p';
  if(l == 2)
    a = 'd';
  if(l == 3)
    a = 'd';
  if(l == 4)
    a = 'f';

  if(n == 1)
    b = '1';
  if(n == 2)
    b = '2';
  if(n == 3)
    b = '3';
  if(n == 4)
    b = '4';
  if(n == 5)
    b = '5';
  //cout << b << a;
 }
 cout << endl;

    
//Termine cinetico
 double sum_cin=0;
 for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k){
  //Non necessario ma velocizza
  if (d.psi[k].Get_n_pieni()>0){     
  //Somme su i e j
   for (unsigned int i=0; i<d.size(); ++i)
    for (unsigned int j=0; j<d.size(); ++j)
     sum_cin+=d.psi[k].Get_n_pieni()*T[ind(i,j)]*
            	d.psi[k].Y_rad_D[i]*d.psi[k].Y_rad_D[j];  
   }
  }
 cout << "Termine cinetico = "<< sum_cin << endl;
 out << sum_cin << "     ";
 E_var+= sum_cin;
 sum_autovalori+= sum_cin;
	

//Termine centrifugo
 unsigned int l;
 double sum_centr=0;

 for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k){      
  //Non necessario ma velocizza
  if (d.psi[k].Get_n_pieni()>0){      
   //Valore di l
   l= d.psi[k].Get_l();   	   
   //cout << "Autofunzione " <<k << " l " << l << endl;     
   //Somme su i
   for (unsigned int i=0; i<d.size(); ++i)
    sum_centr+= d.psi[k].Get_n_pieni()*l*(l+1)/(2*pow(d.X[i],2))
            	*d.psi[k].Y_rad_D[i]*d.psi[k].Y_rad_D[i];
       }
     }

    cout << "Termine centrifugo = "<< sum_centr << endl;
    out << sum_centr << "     ";
    E_var+=sum_centr;
	sum_autovalori+= sum_centr;
	


//Termine potenziale esterno.
  double sum_nucleo=0;
  for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k){     
      //Non necessario ma velocizza
      if (d.psi[k].Get_n_pieni()>0){
      
        //Somme su i
        for (unsigned int i=0; i<d.size(); ++i)
        	sum_nucleo+= d.psi[k].Get_n_pieni()*V_0[i]*d.psi[k].Y_rad_D[i]*d.psi[k].Y_rad_D[i];
       }
     }

  cout << "Termine nucleo = "<< sum_nucleo << endl;
    out << sum_nucleo << "     ";
    E_var+=sum_nucleo;
	sum_autovalori+= sum_nucleo;


//Termine diretto
 double sum_diretto=0;
 for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k){      
      //Non necessario ma velocizza
      if (d.psi[k].Get_n_pieni()>0){
        //Somme su i
        for (unsigned int i=0; i<d.size(); ++i)
        	sum_diretto+= 0.5*d.psi[k].Get_n_pieni()*V_diretto.Y[i]*d.psi[k].Y_rad_D[i]*d.psi[k].Y_rad_D[i];   // *0.5
       }
     }

cout << "Termine diretto = "<< sum_diretto << endl;
 out << sum_diretto << "     ";
 E_var+=sum_diretto;
 sum_autovalori+= 2*sum_diretto;  
 

//Termine potenziale scambio.
 double sum_scambio=0.;

 for (unsigned int i=0; i<d.size(); ++i){
  sum_scambio+= D[i]*pow(d.Y[i]/(4*M_PI*pow(d.X[i],2.)),4./3.)*pow(d.X[i],2.);   
 }
  sum_scambio = -sum_scambio*(3./4.)*pow(3./M_PI,1./3.)*4*M_PI;
  cout << "Termine scambio = "<< sum_scambio << endl;
  out << sum_scambio << "     ";
  E_var+=sum_scambio;
  sum_autovalori+=sum_scambio;

//Termine potenziale di correlazione.

 double p=1.00;
 double A=0.031091;
 double alpha1=0.21370;
 double beta1=7.5957;
 double beta2=3.5876;
 double beta3=1.6382;
 double beta4=0.49294; 

 Vec_DP app(0.,d.size());
 Vec_DP app2(0.,d.size());
 double rs=0.; 

 double sum_corr=0.;

 for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k){      
  //Non necessario ma velocizza
  if (d.psi[k].Get_n_pieni()>0){
   //Somme su i
   for (unsigned int i=0; i<d.size(); ++i){
     double appo=d.Y[i]/(4.*M_PI*pow(d.X[i],2.));
     if(appo != 0){
      rs=pow(3./(4.*M_PI*appo),1./3.); //Se le densità è nulla non è definito.
      double appo1=1.+alpha1*rs;
      app[i]=2.*A*(beta1*sqrt(rs)+beta2*rs+beta3*pow(rs,3./2.)+beta4*pow(rs,p+1.));
      double appo2=1.+1./app[i];
      sum_corr+= d.psi[k].Get_n_pieni()*(-2.*A*appo1*log(appo2))*pow(d.psi[k].Y_rad_D[i],2.);   
     }
    }
   }
  }
  cout << "Termine correlazione = "<< sum_corr << endl;
  out << sum_corr << "     ";
  E_var+=sum_corr;
  sum_autovalori+=sum_corr;


 out << E_var << "          ";
 out << sum_autovalori<< "  ";
    

//Somma degli autovalori
 sum_autovalori=0;
 for (unsigned int i=0; i<d.Get_n_autofunzioni(); ++i) 
  sum_autovalori+= d.psi[i].Get_n_pieni()*d.psi[i].EigenValue();
  //cout << "Somma autovalori = "<< sum_autovalori << endl;
  out << sum_autovalori << "       ";
  
  for (unsigned int k=0; k<d.Get_n_autofunzioni(); ++k) 
   if (d.psi[k].Get_n_pieni()>0){
    out << d.psi[k].Get_n() <<"  " << d.psi[k].Get_l()<< "  "<< 
    d.psi[k].Get_spin()-0.5 << "  " <<d.psi[k].Get_n_pieni() << "        " ;      
    }
  
  if(condition == 1){ 
    ofstream fileout;
    fileout.open("Energia_variazionali_finali.dat",ios::app);
    fileout << Z << " " << N <<" " << E_var << " " << sum_cin << " " << sum_centr << " " << sum_nucleo << " " << sum_diretto << " " << sum_scambio << " " << sum_corr << endl;
    fileout.close();
   }
        
    return E_var;

}
