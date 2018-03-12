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
TEST_CASE("Test thomasAlgorithm_diff", "[ODESolver]"){
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
    auto result=odesolver::thomasAlgorithm_diff(lower, main, upper, solution);
    std::array<double, 3> expected;
    expected[0]=-1.5714286;
    expected[1]=.8571429;
    expected[2]=-1.0714286;
    for(int i=0; i<result.size();++i){
        REQUIRE(result[i]==Approx(expected[i]));
    }
} 
TEST_CASE("Test solveODE", "[ODESolver]"){
    auto fn2=[](const auto& x){
        return 1.5;
    };
    auto fn1=[](const auto& x){
        return 5.0;
    };
    auto fn=[](const auto& x){
        return 1.5;
    };
    double initCondLower=0;
    double initCondUpper=1;
    double xMin=0;
    double xMax=1;
    constexpr int N=100;
    double dx=(xMax-xMin)/(double)(N+1);
    auto result=odesolver::solveODE(fn2, fn1, fn, initCondLower, initCondUpper, xMin, xMax, N);
    auto expectedFunction=[&](const auto& x){
        double c2=1.0/(exp(-1.0/3.0)-exp(-3.0));
        double c1=-c2;
        return c1*exp(-3.0*x)+c2*exp(-(1.0/3.0)*x);
    };
    for(int i=0; i<result.size();++i){
        //std::cout<<result[i]<<", "<<expectedFunction(xMin+(i+1)*dx)<<std::endl;
        REQUIRE(result[i]==Approx(expectedFunction(xMin+(i+1)*dx)).epsilon(.0001));
    }
} 
TEST_CASE("Test solveODE_diff", "[ODESolver]"){
    auto fn2=[](const auto& x){
        return 1.5;
    };
    auto fn1=[](const auto& x){
        return 5.0;
    };
    auto fn=[](const auto& x){
        return 1.5;
    };
    double initCondLower=0;
    double initCondUpper=1;
    double xMin=0;
    double xMax=1;
    constexpr int N=100;
    double dx=(xMax-xMin)/(double)(N+1);
    auto result=odesolver::solveODE_diff(fn2, fn1, fn, initCondLower, initCondUpper, xMin, xMax, N);
    auto expectedFunction=[&](const auto& x){
        double c2=1.0/(exp(-1.0/3.0)-exp(-3.0));
        double c1=-c2;
        return c1*exp(-3.0*x)+c2*exp(-(1.0/3.0)*x);
    };
    for(int i=0; i<result.size();++i){
        //std::cout<<result[i]<<", "<<expectedFunction(xMin+(i+1)*dx)<<std::endl;
        REQUIRE(result[i]==Approx(expectedFunction(xMin+(i+1)*dx)).epsilon(.0001));
    }
} 
