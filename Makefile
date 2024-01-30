ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --glibs) -lRooFit -lRooFitCore -lMinuit

all: Fake

Fake: plotFakenew.cc
	g++ -Wall -Wextra -O3 -o $@ plotFakenew.cc $(ROOTFLAGS) $(ROOTLIBS) 
FakeLocal: plotFakeLocal.cc
	g++ -Wall -Wextra -O3 -o $@ plotFakeLocal.cc $(ROOTFLAGS) $(ROOTLIBS) 
clean:
	rm -f Fake *~

