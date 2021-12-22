CC	      = gcc
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /usr/local/include/gsl/
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lm -lstdc++ 
BAS		= basicmodel.o io.o utilities.o distributions.o likelihoods.o diagnostics.o

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o run_basic  $(LD_FLAGS)
basicmodel.o: basicmodel.cpp
	$(CC) $(CC_FLAGS) -c basicmodel.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp
distributions.o: distributions.cpp
	 $(CC) $(CC_FLAGS) -c distributions.cpp
likelihoods.o: likelihoods.cpp
	$(CC) $(CC_FLAGS) -c likelihoods.cpp
diagnostics.o: diagnostics.cpp
	 $(CC) $(CC_FLAGS) -c diagnostics.cpp
