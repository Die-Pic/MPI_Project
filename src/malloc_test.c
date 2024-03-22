#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

int *test = NULL;

int main(int argc, char const *argv[])
{
    test = malloc(sizeof(int)*0);
    free(test);

    return 0;
}
