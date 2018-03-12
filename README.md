| [Linux][lin-link] | [Windows][win-link] | [Codecov][cov-link] |
| :---------------: | :-----------------: | :-------------------: |
| ![lin-badge]      | ![win-badge]        | ![cov-badge]          |

[lin-badge]: https://travis-ci.org/phillyfan1138/ODESolver.svg?branch=master "Travis build status"
[lin-link]:  https://travis-ci.org/phillyfan1138/ODESolver "Travis build status"
[win-badge]: https://ci.appveyor.com/api/projects/status/g1ii4a3trn7n40c2?svg=true
 "AppVeyor build status"
[win-link]:  https://ci.appveyor.com/project/phillyfan1138/odesolver "AppVeyor build status"
[cov-badge]: https://codecov.io/gh/phillyfan1138/ODESolver/branch/master/graph/badge.svg
[cov-link]:  https://codecov.io/gh/phillyfan1138/ODESolver


# ODESolver
This repository can solve second order ODEs: h(x)f''(x)+g(x)f'(x)+c(x)f(x)=0.  The user must specify 
* the functions h(x), g(x), and c(x)
* two boundary conditions (at f(xmin) and f(xmax)).
* the domain of x (xmin and xmax)
* the number of discrete points to evaluate the function at.  