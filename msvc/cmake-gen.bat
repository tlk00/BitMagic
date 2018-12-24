REM Helper file to generate MSVC build environment with cmake

REM Assign optimization options BMSSE42OPT or BMAVX2OPT
REM cmake -DBMOPTFLAGS:STRING=none ..
REM cmake -DBMOPTFLAGS:STRING=BMSSE42OPT ..
REM cmake -DBMOPTFLAGS:STRING=BMAVX2OPT ..

REM Generate Visual Studio project files
cmake .. -G "Visual Studio 15 2017 Win64"
