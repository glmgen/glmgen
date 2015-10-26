

ROOT=../../
cd $ROOT && \
sudo R CMD INSTALL glmgen && cd $ROOT/glmgen && \
cp src/glmgen.so src/libglmgen.so && \
gcc -g tests/pois_test.c -I inst/include/ -Lsrc/ -lm -lglmgen
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOT/src
