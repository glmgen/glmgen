CC=gcc-4.9 
CFLAGS=-O3 -Wall -Wextra -ansi -std=c89 -pedantic
CFLAGS2=-O3
OBJ=obj/*.o
IDIR=../include/

all:
	mkdir -p lib
	mkdir -p obj
	cd obj; ${CC} ${CLFAGS2} -c -fPIC ../src/csparse/*.c -I${IDIR}
	cd obj; ${CC} ${CFLAGS}  -c -fPIC ../src/utils/*.c -I${IDIR}
	cd obj; ${CC} ${CFLAGS}  -c -fPIC ../src/tf/*.c -I${IDIR}
	${CC} -shared -o lib/libglmgen.so ${OBJ}

docs:
	doxygen Doxyfile

ex:
	${CC} ${CFLAGS} -c -fPIC ./examples/test_utils.c -I./examples -I./src/utils -I./src/tf -I./include/ -o examples/test_utils.o

clean:
	rm obj/*.o examples/*.o

cleanall:
	rm obj/*.o
	rm lib/*.so
	rm -rf examples/bin/*
