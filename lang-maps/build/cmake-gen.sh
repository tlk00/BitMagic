#/bin/sh
# Generate out of the make build using cmake
#
# The build system supports
#  1. default build
#  2. BMSSE42OPT  - SSE4.2 optimizations
#  3. BMAVX2OPT   - AVX2 optimizations
#
# Specific build can be activated by setting cmake cache variable
# Example: cmake -DBMOPTFLAGS:STRING=BMSSE42OPT


#cmake -DBMOPTFLAGS:STRING=BMSSE42OPT ..
#cmake -DBMOPTFLAGS:STRING=BMAVX2OPT ..
cmake ..
