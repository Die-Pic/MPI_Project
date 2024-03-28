MATH_FLAG = -lm

hello: src/fish.c
	@mpicc src/fish.c -o hello $(MATH_FLAG)

test: src/fish_test.c
	@mpicc src/fish.c -i hello $(MATH_FLAG)

clean:
	@rm hello