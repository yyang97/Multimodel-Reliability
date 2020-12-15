
CC = g++

ODIR = ../obj
OPT = -std=c++11 -fopenmp -Ofast

OBJ = model.o utils.o brent.o 
LIBS = -llapack
BOOST = -Iboost

all: ../DOLmod_MCMC ../CanADM ../USADM ../CanADM_snow_Van ../CanADM_wind ../CanADM_snow_Queb ../USADM_snow_Van ../USADM_snow_Queb ../USADM_wind


%.o: %.cpp %.hpp
	$(CC) -c $(OPT) -o $@ $< $(BOOST) 

mcmc_abc.o: mcmc_abc.cpp
	$(CC) -c $(OPT) -o $@ $< $(BOOST)

../DOLmod_MCMC: mcmc_abc.o $(OBJ)
	$(CC) -g $(OPT) -o $@ $< $(OBJ) $(LIBS) $(BOOST)

../CanADM: solveCanODE.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

../CanADM_snow_Van: solveCanODE_snow_Van.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

../CanADM_snow_Queb: solveCanODE_snow_Queb.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

../CanADM_wind: solveCanODE_wind.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

../USADM: solveUSODE.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

../USADM_snow_Van: solveUSODE_snow_Van.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)
../USADM_snow_Queb: solveUSODE_snow_Queb.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)
../USADM_wind: solveUSODE_wind.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

clean:
	rm -f  *.o *~ 


