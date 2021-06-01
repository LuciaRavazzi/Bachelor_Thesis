#include "funzioni.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>


using namespace std;


Dati::Dati() {
    
    ifstream file("input.dat");
    string s;
    if (!file) cerr << "Unable to find file "  << endl;
    else
    {
    file >> N_log;
    file >> Rc; 
    cout << "Raggio cut-off " << Rc << endl;
    file >> Z; 
    cout << "Numero atomico " << Z << endl;
    file >> N;        
    cout << "Numero elettroni " << N << endl;
    file >> err;
    file >> controllo_cicli;
    file >> alpha;                 
    }
    file.close();
}


//genero la mesh
//nota: all'interno della funzione chiamo Rc "max"
Vec_DP mesh(int N_log, double Max, unsigned int Z) {
    
    //calcolo min
    double Min = 0.01/Z; // 1E-6?
    
    //calcolo p della mesh logaritmica
    double p = pow( (Max/Min),(1./(double)(N_log)) ) -1;
    cout << "p = " << p << endl;
    
    //calcolo il numero di punti della mesh lineare tra 0 e min
    //posso usare min(1./p-1,10.) o 1./p-1
    unsigned int N_lin = 1./p-1;  //parte intera
    cout << "I punti della mesh lineare sono " << N_lin << endl <<endl;
    
    //numero di punti della mesh totale
    unsigned int N_tot = N_lin+N_log;
    
    
    vector<double> vett;	
    
    //passo mesh lineare
    double p_lineare = Min/(double)(N_lin+1);
    
    //inserisco mesh lineare
    double val=p_lineare; 	    
    for(int i=0; i< N_lin; ++i){
        vett.push_back(val); 
        val+= p_lineare;
    }
    
    
    //inserisco mesh logaritmica
    val=Min; // posizione di min in vett è N_lin 
    for(int i= N_lin; i<N_tot; ++i){
        vett.push_back(val); 
        val=(1+p)*vett[i];    
    }
       
        
    Vec_DP ris(vett.size());
    
    for (unsigned int i=0; i<N_tot; ++i){
        ris[i]=vett[i];
    }
        
    return ris;    
}


Vec_DP mesh_lin(int N, double Max, unsigned int Z) {
    
    //calcolo min
    double Min = 0.01/Z;
    
    //calcolo il passo p della mesh lineare
    double p = (Max-Min)/(double)N;
    cout << "Il passo è " << p << endl;
    
    vector<double> vett;
    
    //inserisco mesh lineare
    for(int i=0; i<N; ++i){
    	double val=Min+i*p; 
        vett.push_back(val); 
    }
      
    Vec_DP ris(vett.size());
    for (unsigned int i=0; i<N; ++i){
        ris[i]=vett[i];
    }
    
    
    //Controllo
    cout <<"Numero di punti:  " << ris.size() << endl;
    for (unsigned int i=0; i<ris.size(); ++i) 
      cout << i << "      "  <<ris[i]   << endl;    
    
    return ris;    
}


Vec_DP mesh_log(int N, double Max, unsigned int Z) {
    
    //calcolo min
    double Min = 1E-6;
    
    //calcolo p della mesh logaritmica
    double p = pow( (Max/Min),(1./(double)(N)) ) -1;
    cout << "p= " << p << endl;
    
    vector<double> vett;	
    
    //inserisco mesh logaritmica
    double val=Min;  
    for(int i=0; i<N; ++i){
        vett.push_back(val); 
        val=(1+p)*vett[i];    
    }
       
        
    Vec_DP ris(vett.size());
    
    for (unsigned int i=0; i<N; ++i){
        ris[i]=vett[i];
    }
        
    
    //Controllo
    cout <<"Numero di punti:  " << ris.size() << endl;
    for (unsigned int i=0; i<ris.size(); ++i) 
      cout << i << "      "  <<ris[i]   << endl;    
	
    
    return ris;    
}


double factorial(unsigned int k){
    if(k==0)
        return 1;
    else
        return k*factorial(k-1);
}



void set_factorial(double *fact){
    for(int i=0; i<=30; ++i)
        fact[i]=factorial(i);
}


void create_folder(char *nomefile){
    DIR* dir = opendir(nomefile);
	if (dir){  /* Directory exists. */
	  closedir(dir);
	}
	else if (ENOENT == errno){ /* Directory does not exist. */
	  int status;
	  status = mkdir(nomefile, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); 
	  if(status)
	    cerr<<"Error about creating "<< nomefile <<endl; 
	}
	else{ /* opendir() failed for some other reason. */
	  cerr<<"Error about checking "<< nomefile <<endl;
	}
}





