INCLUDES=  -I../FunctionalUtilities
GCCVAL=g++

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	GCCVAL=g++-7
endif
test:test.o 
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage -g  test.o $(INCLUDES) -o test -fopenmp

test.o: test.cpp ODESolver.h 
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage  -g  -c test.cpp   $(INCLUDES) -fopenmp

clean:
	-rm *.o test *.out *.gcda *.gcno



