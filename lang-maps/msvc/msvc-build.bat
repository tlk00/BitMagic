REM MSVC Win64 buld script (requires cmake)
REM

set msdevenv="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\IDE\devenv.com"

mkdir artefact

cmake -DBMOPTFLAGS:STRING=none ..
%msdevenv% /clean Release bmlangmap.sln
%msdevenv% /build Release /out build.log bmlangmap.sln
copy /B /V /Y ..\build\bin\Release\bmcpuid.dll .\artefact\
copy /B /V /Y ..\build\bin\Release\bm-dll.dll .\artefact\

cmake -DBMOPTFLAGS:STRING=BMSSE42OPT ..
%msdevenv% /clean Release bmlangmap.sln
%msdevenv% /build Release /out build.log bmlangmap.sln

copy /B /V /Y ..\build\bin\Release\bm-dll-sse42.dll .\artefact\

cmake -DBMOPTFLAGS:STRING=BMAVX2OPT ..
%msdevenv% /clean Release bmlangmap.sln
%msdevenv% /build Release /out build.log bmlangmap.sln

copy /B /V /Y ..\build\bin\Release\bm-dll-avx2.dll .\artefact\