To build BitMagic C++ using make:

cmake ..
make

A few other variants to configure build for platform specific optimizations:

cmake -DBMOPTFLAGS:STRING=BMSSE2OPT ..
cmake -DBMOPTFLAGS:STRING=BMSSE42OPT ..
cmake -DBMOPTFLAGS:STRING=BMAVX2OPT ..
cmake -DBMOPTFLAGS:STRING=none ..