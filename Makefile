R_PKG_VERSION=0.0.2
PREFIX=

R_DIR=R_pkg
C_DIR=c_lib/glmgen

all:
	git checkout gh-pages
	git rebase master
	cd ${R_DIR}; ${PREFIX} R CMD build glmgen
	cd ${R_DIR}; ${PREFIX} R CMD CHECK --as-cran glmgen_${R_PKG_VERSION}.tar.gz
	cd ${R_DIR}; ${PREFIX} mv glmgen.Rcheck/glmgen-manual.pdf ../..
	cd ${R_DIR}; ${PREFIX} rm -rf glmgen.Rcheck
	cd ${R_DIR}; ${PREFIX} rm glmgen_${R_PKG_VERSION}.tar.gz

	cd ${C_DIR}; doxygen Doxyfile
	cd ${C_DIR}; mv html ..
	cd ${C_DIR}; rm -rf latex

	git add *
	git commit -m "auto commit documentation"
	git push origin gh-pages
	git checkout master
