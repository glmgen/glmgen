

ROOT=../../
cd $ROOT && \
sudo R CMD INSTALL glmgen && cd glmgen && \
cp src/glmgen.so src/libglmgen.so && cd tests && \
gcc -g ctest.c -I../inst/include/ -L../src/ -lm -lglmgen
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOT/src
