#/bin/sh
# Generate out of the source xcode project using cmake
#cmake -DBMOPTFLAGS:STRING=BMSSE2OPT ..
#cmake -DBMOPTFLAGS:STRING=BMSSE42OPT ..
#cmake -DBMOPTFLAGS:STRING=BMAVX2OPT ..
cmake -DBMOPTFLAGS:STRING=none ..
cmake .. -GXcode
