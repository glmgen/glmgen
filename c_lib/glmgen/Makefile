all:
	cd obj; gcc -ansi -Wall -c -fPIC ../src/csparse/*.c -I../src/csparse
	cd obj; gcc -ansi -Wall -c -fPIC ../src/utils/*.c -I../src/csparse -I../src/utils
	cd obj; gcc -ansi -Wall -c -fPIC ../src/tf/*.c -I../src/csparse -I../src/utils -I../src/tf
	gcc -shared -o lib/libglmgen.so obj/*.o

ex:
	gcc -ansi -Wall -I./src/utils -I./src/tf -I./include/ -L./lib -o examples/bin/test examples/test.c -lglmgen -lm
	gcc -ansi -Wall -I./src/utils -I./src/tf -I./include/ -L./lib -o examples/bin/test_mult examples/test_multipliers.c -lglmgen -lm

clean:
	rm obj/*.o

cleanall:
	rm obj/*.o
	rm lib/*.so
	rm examples/bin/*
