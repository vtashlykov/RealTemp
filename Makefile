#
CXX = g++
#CXX = gcc
CXXFLAGS = -Wextra -g3 -O3 -std=c++14 -fopenmp
OUT	 = RealTemp
#
RealTemp: RealTemp.o RadarData.o RadarTime.o InverseProblemsolver.o
	$(CXX) $(CXXFLAGS) -o $(OUT) RealTemp.o RadarData.o RadarTime.o InverseProblemsolver.o -lgsl -lgslcblas -lnlopt -lm
clean:
	rm *.o $(OUT)
