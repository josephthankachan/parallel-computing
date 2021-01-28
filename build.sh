#!/bin/bash

set -e
# echo '###GCC version info###'
# gcc --version
# echo '###GCC version info end###'
ESDK=${EPIPHANY_HOME}
ELIBS="-L ${ESDK}/tools/host/lib"
EINCS="-I ${ESDK}/tools/host/include"
ELDF=${ESDK}/bsps/current/fast.ldf

SCRIPT=$(readlink -f "$0")
EXEPATH=$(dirname "$SCRIPT")
cd $EXEPATH

# Create the binaries directory
mkdir -p bin/

# if [ -z "${CROSS_COMPILE+xxx}" ]; then
# case $(uname -p) in
# 	arm*)
# 		# Use native arm compiler (no cross prefix)
# 		CROSS_COMPILE=
# 		;;
# 	   *)
# 		# Use cross compiler
# 		CROSS_COMPILE="arm-linux-gnueabihf-"
# 		;;
# esac
# fi

# Build all 
gcc -O2 -std=c99 src/main.c src/Sequence_main.c src/Give_Qr_matri.c -o bin/main ${EINCS} ${ELIBS} -le-hal -le-loader -lpthread -lm
e-gcc -Os -std=c99 -falign-loops=8 -falign-functions=8 -Wall  -T ${ELDF} src/epi_seq.c -o bin/epi_seq.elf -le-lib -lm -g
e-gcc -Os -std=c99 -falign-loops=8 -falign-functions=8 -Wall  -T ${ELDF} src/epi_two_core.c -o bin/epi_two_core.elf -le-lib -lm -g
e-gcc -Os -std=c99 -falign-loops=8 -falign-functions=8 -Wall  -T ${ELDF} src/epiphany_four_core_util.c -o bin/epiphany_four_core_util.elf -le-lib -lm -g
