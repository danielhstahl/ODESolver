set mypath=%cd%
cmd /k "cd \Program Files (x86)\Microsoft Visual Studio 14.0\VC & vcvarsall amd64 & cd /d %mypath% & cl /EHsc  /O2 test.cpp /I ..\FunctionalUtilities  & exit 0"