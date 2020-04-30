################################################################################
# Libraries and parse arguments
################################################################################
CC = g++
DEBUG =
OPTIMIZE = -mcmodel=medium -lm -g -Wno-deprecated -Wno-write-strings
LIBS = -lm -fopenmp -lgomp
LIBS = -lfftw3
LIBS_GSL = -lgsl -lgslcblas
#APP_STL = stlport_static
################################################################################

################################################################################
# Compiling and Dependencies
################################################################################

redpower: power_spectrum.o
	$(CC) power_spectrum.o -o redpower $(LIBS) $(DEBUG) $(OPTIMIZE)

power_spectrum.o: power_spectrum.cpp
	$(CC) -c power_spectrum.cpp -g $(LIBS) $(OPTIMIZE)
	
clean:
	rm *.o redpower

