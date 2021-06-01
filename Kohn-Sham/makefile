
CXX = g++

CXXFLAGS = $(PRE_FLAG) -O3 -ansi -funroll-loops -fexpensive-optimizations  #-static
#CXXFLAGS = $(PRE_FLAG) -ansi -pedantic -Wall -Winline -DCHECK_LIMITS -g -ggdb # move below to use


LDFLAGS = -llapack 

hartree-fock: main.o lapack.o funzioni.o my_nr.o spline.o splint.o function.o wave_function.o density.o potential.o e_var.o evaluate.o V_Diretto.o 
	$(CXX) $(CXXFLAGS) main.o lapack.o funzioni.o my_nr.o spline.o splint.o function.o wave_function.o density.o potential.o e_var.o evaluate.o V_Diretto.o -o hartree-fock $(LDFLAGS)


main.o: lapack.h funzioni.h density.h nr.h potential.h e_var.h evaluate.h main.cpp V_Diretto.h 
	$(CXX) $(CXXFLAGS) -c main.cpp 

lapack.o: nr.h lapack.h lapack.cpp
	$(CXX) $(CXXFLAGS) -c lapack.cpp 

funzioni.o: nr.h funzioni.h funzioni.cpp
	$(CXX) $(CXXFLAGS) -c funzioni.cpp

my_nr.o: nr.h my_nr.h my_nr.cpp
	$(CXX) $(CXXFLAGS) -c my_nr.cpp

spline.o: nr.h spline.cpp
	$(CXX) $(CXXFLAGS) -c spline.cpp

splint.o: nr.h splint.cpp
	$(CXX) $(CXXFLAGS) -c splint.cpp

function.o: nr.h function.h function.cpp
	$(CXX) $(CXXFLAGS) -c function.cpp

wave_function.o: nr.h function.h wave_function.h wave_function.cpp
	$(CXX) $(CXXFLAGS) -c wave_function.cpp

density.o: nr.h lapack.h funzioni.h function.h wave_function.h my_nr.h density.h density.cpp
	$(CXX) $(CXXFLAGS) -c density.cpp

potential.o: nr.h function.h density.h potential.h potential.cpp
	$(CXX) $(CXXFLAGS) -c potential.cpp

e_var.o: funzioni.h wave_function.h density.h potential.h e_var.h e_var.cpp V_Diretto.h 
	$(CXX) $(CXXFLAGS) -c e_var.cpp

evaluate.o: nr.h lapack.h funzioni.h function.h wave_function.h my_nr.h density.h potential.h evaluate.cpp
	$(CXX) $(CXXFLAGS) -c evaluate.cpp

V_Diretto.o: V_Diretto.cpp V_Diretto.h nr.h function.h density.h
	$(CXX) $(CXXFLAGS) -c V_Diretto.cpp


clean:
	rm hartree-fock *.o;
	clear









