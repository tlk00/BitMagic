#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

echo "DEBUG"
#./perf_debug || exit 1

echo
echo
echo RELEASE

./perf_release || exit 1

echo
echo
echo NEON

./perf_release_neon || exit 1

echo
echo "BitMagic Perf Test DONE"
