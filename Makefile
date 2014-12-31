R_PKG_VERSION=0.0.2
CC=gcc-4.9

R_DIR=R_pkg
C_DIR=c_lib/glmgen
CFLAGS=-O3 -Wall -Wextra -ansi -std=c89 -pedantic
CFLAGS2=-O3
OBJ=obj/*.o
IDIR=../include/

all:
	cd ${R_DIR}; R CMD build glmgen
	cd ${R_DIR}; R CMD INSTALL glmgen_${R_PKG_VERSION}.tar.gz
	cd ${R_DIR}; R CMD CHECK --as-cran glmgen_${R_PKG_VERSION}.tar.gz
	cd ${R_DIR}; rm -rf glmgen.Rcheck
	cd ${R_DIR}; rm glmgen_${R_PKG_VERSION}.tar.gz

	cp R_pkg/glmgen/src/tf/*.c c_lib/glmgen/src/tf/
	cp R_pkg/glmgen/src/utils/*.c c_lib/glmgen/src/utils/
	cp R_pkg/glmgen/inst/include/*.h c_lib/glmgen/include

	cd ${C_DIR}; mkdir -p lib
	cd ${C_DIR}; mkdir -p obj
	cd ${C_DIR}/obj; ${CC} ${CLFAGS2} -c -fPIC ../src/csparse/*.c -I${IDIR}
	cd ${C_DIR}/obj; ${CC} ${CFLAGS}  -c -fPIC ../src/utils/*.c -I${IDIR}
	cd ${C_DIR}/obj; ${CC} ${CFLAGS}  -c -fPIC ../src/tf/*.c -I${IDIR}
	cd ${C_DIR}; ${CC} -shared -o lib/libglmgen.so ${OBJ}