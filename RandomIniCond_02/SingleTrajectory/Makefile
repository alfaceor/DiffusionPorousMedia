programs=rw3D_RegularPorous__time_1e8__RandomIniCond_02 rw3D_RegularPorous__time_1e7__RandomIniCond_02 rw3D_RegularPorous__time_1e6__RandomIniCond_02 rw3D_RegularPorous__time_1e4__RandomIniCond_02 rw3D_Histo__RandomIniCond_02  parallel__time_1e4__openmp parallel__time_1e7__openmp parallel__time_1e8__openmp parallel__time_1e9__openmp rw3D_Geometry
all:$(programs)

LIBS= -lgsl -lgslcblas -lm -lboost_program_options

rw3D_Geometry:rw3D_Geometry.cpp
	g++ -O3 -o rw3D_Geometry rw3D_Geometry.cpp $(LIBS)

parallel__time_1e4__openmp:rw3D_RegularPorous__time_1e4__RandomIniCond_02_openmp.cpp
	g++ -O3 -o parallel__time_1e4__openmp -fopenmp  rw3D_RegularPorous__time_1e4__RandomIniCond_02_openmp.cpp $(LIBS)

parallel__time_1e7__openmp:rw3D_RegularPorous__time_1e7__RandomIniCond_02_openmp.cpp
	g++ -O3 -o parallel__time_1e7__openmp -fopenmp  rw3D_RegularPorous__time_1e7__RandomIniCond_02_openmp.cpp $(LIBS)

parallel__time_1e8__openmp:rw3D_RegularPorous__time_1e8__RandomIniCond_02_openmp.cpp
	g++ -O3 -o parallel__time_1e8__openmp -fopenmp  rw3D_RegularPorous__time_1e8__RandomIniCond_02_openmp.cpp $(LIBS)

parallel__time_1e9__openmp:rw3D_RegularPorous__time_1e9__RandomIniCond_02_openmp.cpp
	g++ -O3 -o parallel__time_1e9__openmp -fopenmp  rw3D_RegularPorous__time_1e9__RandomIniCond_02_openmp.cpp $(LIBS)

rw3D_RegularPorous__time_1e8__RandomIniCond_02:rw3D_RegularPorous__time_1e8__RandomIniCond_02.cpp
	g++ -O3 -o rw3D_RegularPorous__time_1e8__RandomIniCond_02 rw3D_RegularPorous__time_1e8__RandomIniCond_02.cpp $(LIBS)

rw3D_RegularPorous__time_1e7__RandomIniCond_02:rw3D_RegularPorous__time_1e7__RandomIniCond_02.cpp
	g++ -O3 -o rw3D_RegularPorous__time_1e7__RandomIniCond_02 rw3D_RegularPorous__time_1e7__RandomIniCond_02.cpp $(LIBS)

rw3D_RegularPorous__time_1e6__RandomIniCond_02:rw3D_RegularPorous__time_1e6__RandomIniCond_02.cpp
	g++ -O3 -o rw3D_RegularPorous__time_1e6__RandomIniCond_02 rw3D_RegularPorous__time_1e6__RandomIniCond_02.cpp $(LIBS)

rw3D_RegularPorous__time_1e4__RandomIniCond_02:rw3D_RegularPorous__time_1e4__RandomIniCond_02.cpp
	g++ -O3 -o rw3D_RegularPorous__time_1e4__RandomIniCond_02 rw3D_RegularPorous__time_1e4__RandomIniCond_02.cpp $(LIBS)

rw3D_Histo__RandomIniCond_02:rw3D_Histo__RandomIniCond_02.cpp
	g++ -O3 -o rw3D_Histo__RandomIniCond_02 rw3D_Histo__RandomIniCond_02.cpp $(LIBS)
clean:
	rm  $(programs)
