#include <iostream>
#include <cmath>
#include <fstream>
#include "integrale.h" 
#include "funzioneBase.h" 
#include "integranda.h"
#include "funzioni.h"
#include "main.h"
#include <vector>
#include "nr.h"
#include <iomanip>

using namespace std;


int main(){

    ofstream out;
  if(dati.boolean == 1){ 
    out.open("N(mu)Z="+to_string(dati.Z)+  ".dat");
    double mu=dati.muInf;
    double passo_mu=abs(dati.muInf-dati.muSup)/100;
    cout << "Passo potenziale chimico " << passo_mu << endl;
    do{
      Iniz();
      double a=Nel(mu);
      out << mu << " " << a << endl; 
      cout << mu << " " << a << endl;  
      mu+=passo_mu;
    } while(mu <= dati.muSup);
  
    out.close();
   }

   Iniz();
   cout << "Z " << Z << endl; 
  //Calcolo il numero di particelle con un \mu molto grande negativo. 
  //double NelInf=Nel(dati.muInf);
  //cout << "muInf " << dati.muInf << " NelInf " << NelInf << endl;
  //PrintEnergy();
  //Ripeto ma con \mu molto vicino allo zero (\mu > 0 stati non legati).
  //double NelSup=Nel(dati.muSup);  
  //cout << "muSup " << dati.muSup << " NelSup " << NelSup << endl;
  //PrintEnergy();

  //Implemento il metodo della bisezione.
   for(int N=Z-dati.abbassaN; N<=Z+dati.alzaN; N++){ 
     if (N > 0){
       fileout.open("Density/Density(Z=" + to_string(dati.Z) + ",N=" + to_string(N) + ").dat");
       DP root= NR::rtbis(Nel,dati.muInf,dati.muSup,dati.acc,N);
       cout << " 	*** RADICE ***	" << root << endl;
       //Memorizzo la densità associata al potenziale chimico trovato.
       //Iniz();
       double numerotrovati=Nel(root);
       cout << "numerotrovati " << numerotrovati << endl;
       cout << "Stampo il file" << endl;
       for(int i=0; i<dim; i++){
         fileout << X_fissate[i] << " " << nNew[i] << " " << 4*M_PI*nNew[i]*X_fissate[i]*X_fissate[i] <<  endl;
         //cout << X_fissate[i] << " " << nNew[i] << endl;
       }
       fileout.close();
     }
     if(Z == N) PrintEnergy();
   }


//Nella tesi bisogna scrivere sia dell'oscillazione per ioni.
//Correzione dovuto al fatto che un elettrone non interagisce con se stesso.
return 0;
}

//CREO LA PRIMA DENSITÀ DI PROVA.
void Iniz(){
  for(int i=0; i<dim; i++) nOld[i]=exp(-X_fissate[i]);
  for(int i=0; i<dim-1; i++) diff[i]=X_fissate[i+1]-X_fissate[i];
}

//CALCOLO IL NUMERO DI ELETTRONI A FISSATO \MU.
double Nel(double mu){
 double E_old;
 double E_new;
 double eps;
 Iniz(); //Ad ogni ciclo devo dire la prima densità di prova.
 //ofstream fileout;
 //fileout.open("Densità");
 int counter=0;
 int iter=0;
 do{
   counter++;
   double p=Particles(nOld);
   corrFactor=(p-1.)/p; 
   if(corrFactor < 0) corrFactor=0.;  //Quando ho meno di un elettrone, non ho interazione tra elettroni.
     iter++;
     //fileout << endl << "# " << iter << endl;
     sum=0.;
     E_old = E_new; 
     //cout << "E_old " << E_old << endl;
     for(int i=0; i<dim; i++){
	//Devo considerare che un elettrone non interagisce con se stesso.
       double b=A*(mu-potExt(X_fissate[i],Z,rn)-corrFactor*Hart.trapezi(X_fissate,nOld,dim,X_fissate[i]));
       if(b>0){
         nNew[i]=2*C*(pow(b,1.5));
         } else {
           nNew[i]=0;
         }
      }
      eps=Distanza();
      if(counter%10 == 0) cerr << iter << " " << mu << " " << " Distanza " << eps << " Nel " << Particles(nNew) << " corrFactor " << corrFactor << endl; //Commentare quando tutto funziona.
      for(int i=0; i<dim; i++){
       nOld[i]=0.9*nOld[i]+0.1*nNew[i]; //Il mixing va bene anche da quando b<0. Posso cambiare i parametri, influisce sulla velocità del programma.
       //fileout << X_fissate[i] << " " << nOld[i]*X_fissate[i]*X_fissate[i] << " " << nNew[i]*X_fissate[i]*X_fissate[i] << endl;
       //La nuova densità diventa quella vecchia
      }
      //   PrintEnergy();
      E_new = EKin() + EPot() + PotEl(); //energia con la densità nuova.
      //cout << "E_new " << E_new << endl;
    } while(eps > dati.toll); //Metti nel file input.
   // fileout.close();
    err=abs(E_new-E_old);
    cout << "Errore " << err << endl;
return Particles(nNew);
}


//TERMINI ENERGIA DI TF.
double EKin(){ 
  double h=(3./5.)*0.5*pow(C,-2./3.)*4*M_PI;
  double inte=0.;
  double FT=5./3.;
  for(int i=0; i<dim-1; i++) 
    inte += (pow(nOld[i+1],FT)*pow(X_fissate[i+1],2)+pow(nOld[i],FT)*pow(X_fissate[i],2))*diff[i];

  return (inte/2.)*h;
}

double EPot(){
  double inte=0.;
  for(int i=0; i<dim-1; i++) 
   inte += (potExt(X_fissate[i+1],Z,rn)*nOld[i+1]*pow(X_fissate[i+1],2)+potExt(X_fissate[i],Z,rn)*nOld[i]*pow(X_fissate[i],2))*diff[i];

  return 4*M_PI*(inte/2.);
}

double PotEl(){
  double sum=0.;
  double* inte1_g=new double[dim];
  double* inte2_g=new double[dim];
 

//Primo intregrale.
  for(int k=0; k<dim; k++){ 
    sum=0.;
    for(int j=0; j<k; j++) sum += (pow(X_fissate[j+1],2)*nOld[j+1]+pow(X_fissate[j],2)*nOld[j])*diff[j];
    inte1_g[k]=sum/2.;
    //cout << "inte1_g[k] " <<  inte1_g[k] << endl;
  } 

  double inte1=0.;
  for(int j=0; j<dim-1; j++){ 
    inte1 += (X_fissate[j+1]*nOld[j+1]*inte1_g[j+1]+X_fissate[j]*nOld[j]*inte1_g[j])*diff[j]; 
   // cout << "inte1" << inte1 << endl;
   }

//Secondo integrale.
 sum=0.;
 for(int k=0; k<dim; k++){ 
   sum=0.;
    for(int j=k; j<dim-1; j++) sum += (X_fissate[j+1]*nOld[j+1]+X_fissate[j]*nOld[j])*diff[j];
    inte2_g[k]=sum/2.;
     //cout << "inte2_g[k] " << inte2_g[k] << endl;
   }
 
  
   double inte2=0.;
   for(int j=0; j<dim-1; j++){ 
    inte2 += (pow(X_fissate[j+1],2)*nOld[j+1]*inte2_g[j+1]+pow(X_fissate[j],2)*nOld[j]*inte2_g[j])*diff[j]; 
    //cout << "inte2 " << inte2 << endl; 
   }


   return (corrFactor*8*pow(M_PI,2)*(inte1+inte2))/2.;    
}

double potExt(double r, int Z, double rn){ //Nucleo non puntiforme.
  if(r < rn){
    return -(3.*Z)/(2.*rn)+(Z*pow(r,2))/(2*pow(rn,3));
  } else {
    // if( r < dati.rSfera){
      return -Z/r; 
     //} else {
      // return -Z/r+0.001*pow(r-dati.rSfera,2); //Oltre ad un certo r, aggiungo una parabola.
   // }
  }
}

//STAMPO L'ENERGIA DI TF.
void PrintEnergy(){
 ofstream fileout;
 fileout.open("Energy.dat", ios::app);
 double ekin=EKin();
 double ep=EPot();
 double EHa=PotEl();

 cout << "Ekin     " << ekin << endl;
 cout << "PotExt   " << ep << endl;
 cout << "Hartree  " << EHa << endl;
 cout << "ETot     " << ekin+ep+EHa << endl;
 fileout << dati.Z << " " << ekin << " " << ep << " " << EHa << " " << ekin+ep+EHa << " " << err <<  endl;  
}

//METODO DEI TRAPEZI PER IL CALCOLO DEL NUMERO DI PARTICELLE.
double Particles(double* n){ //Calcola il numero di particelle della nuova funzione densità.
  part=0.;
  for(int i=0; i<dim-1; i++) 
    part+=(n[i+1]*pow(X_fissate[i+1],2)+n[i]*pow(X_fissate[i],2))*diff[i];

  return part*((4*M_PI)/2.);
}

//DISTANZA TRA LA VECHHIA E LA NUOVA DENSITÀ.
double Distanza(){
  double sum=0.;

  for(int i=0; i<dim-1; i++) 
    sum+=(pow(nNew[i+1]-nOld[i+1],2)*pow(X_fissate[i+1],2)+pow(nNew[i]-nOld[i],2)*pow(X_fissate[i],2))*diff[i];

  return sum*((4*M_PI)/2.);
}

void Print(int N){
  fileout.open("DensityFin(Z=" + to_string(Z) + "/N=" + to_string(N));
  cout << "Stampo il file" << endl;
  for(int i=0; i<dim; i++)
    fileout << X_fissate[i] << " " << nNew[i] << endl;
  
  cout << "Particelle finali: " << Particles(nNew) << endl;

  fileout.close();
}



