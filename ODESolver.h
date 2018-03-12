#ifndef __ODE_H_INCLUDED__
#define __ODE_H_INCLUDED__

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#include <cmath>
#include "FunctionalUtilities.h"

namespace odesolver {

    

    template<typename lowerArray, typename mainArray, typename upperArray, typename solveArray>
    auto thomasAlgorithm(lowerArray&& lower,mainArray&& main, upperArray&& upper, solveArray&& solve){
        auto lastIndexUpper=upper.size();
        upper=futilities::reduce(upper, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){
                return curr/(main[index]-prev*lower[index-1]);
            }
            else{
                return curr/main[index];
            }
        });
        solve=futilities::reduce(solve, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){ 
                return (curr-lower[index-1]*prev)/(main[index]-lower[index-1]*upper[index-1]);
            }
            else {
                return curr/main[index];
            }
        });
        solve=futilities::reduce_reverse(solve, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){
                return curr-upper[lastIndexUpper-index]*prev;
                
            }
            else{ 
                return curr;
            }
        });
        return std::move(solve);
    }
    /**used if the coefficients of the function have different type than the output.  Less efficient since new arrays are created*/
    template<typename lowerArray, typename mainArray, typename upperArray, typename solveArray>
    auto thomasAlgorithm_diff(const lowerArray& lower, const mainArray& main, const upperArray& upper, const solveArray& solve){
        auto lastIndexUpper=upper.size();
        auto modUpper=futilities::reduce_copy(upper, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){
                return curr/(main[index]-prev*lower[index-1]);
            }
            else{
                return curr/main[index];
            }
        });
       
        auto modSolve=futilities::reduce_copy(solve, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){ 
                return (curr-lower[index-1]*prev)/(main[index]-lower[index-1]*modUpper[index-1]);
            }
            else {
                return curr/main[index];
            }
        });
        return futilities::reduce_reverse(modSolve, [&](const auto& prev, const auto& curr, const auto& index){
            if(index>0){
                return curr-modUpper[lastIndexUpper-index]*prev;
            }
            else{ 
                return curr;
            }
        });
    }


    /**Solves ODEs of the form fn2(x)*f''(x)+fn1(x)*f'(x)+fn*f(x)=0*/
    template<typename CoefSecondDeriv, typename CoefFirstDerivative, typename CoefFunction, typename Number, typename Index>
    auto solveODE(
        const CoefSecondDeriv& fn2, 
        const CoefFirstDerivative& fn1, 
        const CoefFunction& fn,
        const Number& initialConditionLower,
        const Number& initialConditionUpper,
        const Number& xMin,
        const Number& xMax,
        const Index& N
    ){
        auto dx=(xMax-xMin)/(double)(N+1); //actual discrete steps are N+2 (N diaganonal+2 conditions)
        auto dxsq=dx*dx;
        auto dx2=dx*2.0;
        auto getXAtIndex=[&](const auto& index){
            return xMin+dx*index;
        };
        auto getUpperCoef=[&](const auto& index){
            auto x=getXAtIndex(index);
            
            return fn2(x)/dxsq+fn1(x)/dx2;
        };
        auto getLowerCoef=[&](const auto& index){
            auto x=getXAtIndex(index);
            return fn2(x)/dxsq-fn1(x)/dx2;
        };
        auto getMainCoef=[&](const auto& index){
            auto x=getXAtIndex(index);
            return fn(x)-fn2(x)*2.0/dxsq;
        }; 

        auto main=futilities::for_each_parallel(0, N, [&](const auto& index){
            return getMainCoef(index+1);//starts one after xMin...xMin is the boundary
        });
        auto upper=futilities::for_each_parallel(0, N-1, [&](const auto& index){
            return getUpperCoef(index+1);//starts one after the xMin...xMin is the boundary
        });
        auto lower=futilities::for_each_parallel(0, N-1, [&](const auto& index){
            return getLowerCoef(index+2);//starts two after the xMin...the one after xMin is the coefficient on the boundary
        });

       

        auto solution=futilities::for_each_parallel(0, N, [&](const auto& index){
            if(index==0){
                return -initialConditionLower*getLowerCoef(index+1);
            }
            else if(index==N-1){
                
                return -initialConditionUpper*getUpperCoef(index+1);
            }
            else{
                return 0.0;
            }
        });
        return thomasAlgorithm(std::move(lower), std::move(main), std::move(upper), std::move(solution));
    }
    /**Solves ODEs of the form fn2(x)*f''(x)+fn1(x)*f'(x)+fn*f(x)=0*/
    template<typename CoefSecondDeriv, typename CoefFirstDerivative, typename CoefFunction, typename Number, typename Index>
    auto solveODE_diff(
        const CoefSecondDeriv& fn2, 
        const CoefFirstDerivative& fn1, 
        const CoefFunction& fn,
        const Number& initialConditionLower,
        const Number& initialConditionUpper,
        const Number& xMin,
        const Number& xMax,
        const Index& N
    ){
        auto dx=(xMax-xMin)/(double)(N+1); //actual discrete steps are N+2 (N diaganonal+2 conditions)
        auto dxsq=dx*dx;
        auto dx2=dx*2.0;
        auto getXAtIndex=[&](const auto& index){
            return xMin+dx*index;
        };
        auto getUpperCoef=[&](const auto& index){
            auto x=getXAtIndex(index);
            return fn2(x)/dxsq+fn1(x)/dx2;
        };
        auto getLowerCoef=[&](const auto& index){
            auto x=getXAtIndex(index);
            return fn2(x)/dxsq-fn1(x)/dx2;
        };
        auto getMainCoef=[&](const auto& index){
            auto x=getXAtIndex(index);
            return fn(x)-fn2(x)*2.0/dxsq;
        }; 

        auto main=futilities::for_each_parallel(0, N, [&](const auto& index){
            return getMainCoef(index+1);//starts one after xMin...xMin is the boundary
        });
        auto upper=futilities::for_each_parallel(0, N-1, [&](const auto& index){
            return getUpperCoef(index+1);//starts one after the xMin...xMin is the boundary
        });
        auto lower=futilities::for_each_parallel(0, N-1, [&](const auto& index){
            return getLowerCoef(index+2);//starts two after the xMin...the one after xMin is the coefficient on the boundary
        });

        auto solution=futilities::for_each_parallel(0, N, [&](const auto& index){
            if(index==0){
                return -initialConditionLower*getLowerCoef(index+1);
            }
            else if(index==N-1){
                return -initialConditionUpper*getUpperCoef(index+1);
            }
            else{
                return 0.0*main[index];//*main[index];//in case its not a double...
            }
        });
        return thomasAlgorithm_diff(lower, main, upper, solution);
    }
    
}

template<int N, typename CoefSecondDeriv, typename CoefFirstDerivative, typename CoefFunction, typename Number>
auto solveODE(
    const CoefSecondDeriv& fn2, 
    const CoefFirstDerivative& fn1, 
    const CoefFunction& fn,
    const Number& xMin,
    const Number& xMax,
    const Number& initialConditionLower,
    const Number& initialConditionUpper
){
    return solveODE(fn2, fn1, fn, xMin, xMax, initialConditionLower, initialConditionUpper, [](const auto& val){
        return val;
    });
}


#endif