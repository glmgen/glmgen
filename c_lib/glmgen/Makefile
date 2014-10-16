all:
	mkdir -p lib
	mkdir -p obj
	cd obj; gcc -g -ansi -Wall -c -fPIC ../src/csparse/*.c -I../src/csparse
	cd obj; gcc -g -ansi -Wall -c -fPIC ../src/utils/*.c -I../src/csparse -I../src/utils
	cd obj; gcc -g -ansi -Wall -c -fPIC ../src/tf/*.c -I../src/csparse -I../src/utils -I../src/tf
	gcc -shared -o lib/libglmgen.so obj/*.o

ex:
	gcc -g -ansi -Wall -c -fPIC ./examples/test_utils.c -I./examples -I./src/utils -I./src/tf -I./include/ -o examples/test_utils.o
	gcc -g -ansi -Wall -I./src/utils -I./src/tf -I./include/ -I./examples -L./lib -o examples/bin/test examples/test.c obj/*.o examples/*.o -lm
	gcc -g -ansi -Wall -I./src/utils -I./src/tf -I./include/ -I./examples -L./lib -o examples/bin/test_mult examples/test_multipliers.c obj/*.o  examples/*.o -lm
	gcc -g -ansi -Wall -I./src/utils -I./src/tf -I./include/ -I./examples -L./lib -o examples/bin/test_admm examples/test_admm.c obj/*.o  examples/*.o -lm

clean:
	rm obj/*.o examples/*.o

cleanall:
	rm obj/*.o
	rm lib/*.so
	rm examples/bin/*
