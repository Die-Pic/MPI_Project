#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef struct fish_s{
    int32_t x, y, z;
    uint16_t speed;
    uint16_t size;
    uint8_t direction;
}fish;

fish* add_fish(fish *queue, uint32_t *index, uint32_t *length, fish toAdd){
    //If queue is full, increase size
    if((*index) >= (*length)){
        puts("realloc");
        queue = (fish*)realloc(queue, (*length) * 2 * sizeof(fish));
        (*length) = (*length)*2;
    }

    queue[(*index)] = toAdd;
    (*index)++;


    return queue;
}

int main(int argc, char const *argv[])
{
    fish* queue = NULL;
    uint32_t index = 0, length = 0;

    queue = malloc(sizeof(fish)*3);
    length = 3;

    for(int i = 0; i < 4; i++){
        fish f = {1, 1, 2, 3, 4, 5};
        queue = add_fish(queue, &index, &length, f);
    }

    for(int i = 0; i < index; i++){
        printf("fish %d\n", queue[i].size);
    }

    uint32_t test = 0;
    test = (test - 1) % 5;
    printf("Test: %d\n", test);

    return 0;
}
