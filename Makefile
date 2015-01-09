R_PKG_VERSION=0.0.2
PREFIX=

R_DIR=R_pkg
C_DIR=c_lib/glmgen
CFLAGS=-O3 -Wall -Wextra -ansi -std=c89 -pedantic
CFLAGS2=-O3
OBJ=obj/*.o
IDIR=../include/

all:
	cd ${R_DIR}; ${PREFIX} R CMD build glmgen

	cd ${C_DIR}; mkdir -p lib
	cd ${C_DIR}; mkdir -p obj
	cd ${C_DIR}/obj; ${CC} ${CLFAGS2} -c -fPIC ../src/csparse/*.c -I${IDIR}
	cd ${C_DIR}/obj; ${CC} ${CFLAGS}  -c -fPIC ../src/utils/*.c -I${IDIR}
	cd ${C_DIR}/obj; ${CC} ${CFLAGS}  -c -fPIC ../src/tf/*.c -I${IDIR}
	cd ${C_DIR}; ${CC} -shared -o lib/libglmgen.so ${OBJ}

clean:
	cd ${R_DIR}; ${PREFIX} rm -rf glmgen_${R_PKG_VERSION}.tar.gz
	cd ${R_DIR}; ${PREFIX} rm -rf glmgen.Rcheck

	cd ${C_DIR}; rm -rf lib
	cd ${C_DIR}; rm -rf obj

doc:
	git checkout gh-pages
	git rebase master
	cd ${R_DIR}; ${PREFIX} R CMD build glmgen
	cd ${R_DIR}; ${PREFIX} R CMD CHECK --as-cran glmgen_${R_PKG_VERSION}.tar.gz
	cd ${R_DIR}; ${PREFIX} mv glmgen.Rcheck/glmgen-manual.pdf ..
	cd ${R_DIR}; ${PREFIX} rm -rf glmgen.Rcheck
	cd ${R_DIR}; ${PREFIX} rm glmgen_${R_PKG_VERSION}.tar.gz

	cd ${C_DIR}; doxygen Doxyfile
	rm -rf html
	cd ${C_DIR}; mv html ../..
	cd ${C_DIR}; rm -rf latex

	git add -f glmgen-manual.pdf
	git add html
	git commit -m "auto commit documentation"
	git push origin gh-pages
	git checkout master
