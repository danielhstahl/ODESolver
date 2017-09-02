#ifndef __ODE_H_INCLUDED__
#define __ODE_H_INCLUDED__

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#include <cmath>
#include "FunctionalUtilities.h"

namespace odesolver {
    /*template<typename Alpha, typename Sigma, int N>
    auto createMatrix(const Alpha& alpha, const Sigma& sigma, ){
        std::array<std::complex<double>, N> main;
        std::array<double, N-1> upper;
        std::array<double, N-1> lower;

    }*/

    template<typename lowerArray, typename mainArray, typename upperArray, typename solveArray>
    auto thomasAlgorithm(const lowerArray& lower,const mainArray& main, const upperArray& upper, solveArray&& solve){
        auto lastIndexUpper=upper.size();
        auto modUpper=futilities::reduce_copy(upper, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){
                return curr/(main[index]-prev*lower[index-1]);
            }
            else{
                return curr/main[index];
            }
        });

        solve=futilities::reduce(solve, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){
                return (curr-lower[index-1]*prev)/(main[index]-lower[index-1]*modUpper[index-1]);
            }
            else {
                return curr/main[index];
            }
        });
        solve=futilities::reduce_reverse(solve, [&](const auto& prev, const auto& curr, const auto& index){
            std::cout<<"index:"<<index<<std::endl;

            if(index>0){
                return curr-modUpper[lastIndexUpper-index]*prev;
                
            }
            else{ 
                return curr;
            }
        });
        return std::move(solve);
    }

}


#endif