CXX=c++
CXXFLAGS=-O2 -Wall -W  -I.

TARGETS=expevol_region expevol_region_breakpoints expevol_region_breakpoints_Ksel 

all: expevol_region.o expevol_region_breakpoints.o expevol_region_breakpoints_Ksel.o
	$(CXX) $(CXXFLAGS) -o expevol_region expevol_region.o $(LDFLAGS) -lsequence -lz -lgsl -lgslcblas -lboost_iostreams
	$(CXX) $(CXXFLAGS) -o expevol_region_breakpoints expevol_region_breakpoints.o $(LDFLAGS) -lsequence -lz -lgsl -lgslcblas -lboost_iostreams
	$(CXX) $(CXXFLAGS) -o expevol_region_breakpoints_Ksel expevol_region_breakpoints_Ksel.o $(LDFLAGS) -lsequence -lz -lgsl -lgslcblas -lboost_iostreams

clean:
	rm -f *.o $(TARGETS)