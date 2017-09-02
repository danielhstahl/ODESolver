#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include <iostream>
#include "ODESolver.h"
#include <complex>
/**
Compile your application with -g, then you'll have debug symbols in the binary file.
Use gdb to open the gdb console.
Use file and pass it your application's binary file in the console.
Use run and pass in any arguments your application needs to start.
Do something to cause a Segmentation Fault.
Type bt in the gdb console to get a stack trace of the Segmentation Fault.*/
TEST_CASE("Test thomasAlgorithm", "[ODESolver]"){
    std::array<double, 3> main;
    std::array<double, 3> solution;
    std::array<double, 2> lower;
    std::array<double, 2> upper;
    main[0]=.3;
    main[1]=.4;
    main[2]=.2;
    lower[0]=.4;
    lower[1]=.6;
    upper[0]=.9;
    upper[1]=.2; 
    solution[0]=.3;
    solution[1]=-.5;
    solution[2]=.3;
    auto result=odesolver::thomasAlgorithm(lower, main, upper, solution);
    std::array<double, 3> expected;
    expected[0]=-1.5714286;
    expected[1]=.8571429;
    expected[2]=-1.0714286;
    for(int i=0; i<result.size();++i){
        REQUIRE(result[i]==Approx(expected[i]));
    }
} 