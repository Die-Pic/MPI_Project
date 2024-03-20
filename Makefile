MATH_FLAG = -lm

hello: src/fish.c
	@mpicc src/fish.c -o hello $(MATH_FLAG)