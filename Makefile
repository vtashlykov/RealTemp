#
CXX = g++
#CXX = gcc
CXXFLAGS = -Wextra -g3 -O3 -std=c++14 -fopenmp
#
RealTemp: RealTemp.o RadarData.o RadarTime.o InverseProblemsolver.o
	$(CXX) $(CXXFLAGS) -o RealTemp RealTemp.o RadarData.o RadarTime.o InverseProblemsolver.o -lgsl -lgslcblas -lnlopt -lm
clean:
	rm *.o
