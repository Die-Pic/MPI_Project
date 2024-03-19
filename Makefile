MATH_FLAG = -lm

hello:
	@mpicc src/fish.c -o hello $(MATH_FLAG)