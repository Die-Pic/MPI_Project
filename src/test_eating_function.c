#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 100           //Number of fish (parameter)
#define L 100           //Dimension of fishtank (parameter)
#define v 1             //Speed of fish (parameter)
#define d 1             //Max eating distance (parameter)
#define t 100           //Time of each step (parameter)
#define MAX_SIZE 100

typedef struct fish{
    int32_t x, y, z;
    uint16_t speed;
    uint16_t size;
    uint8_t direction;
}fish;

int main(int argc, char const *argv[])
{
    int a[10];
    srand(time(NULL));
    
    for(uint32_t i = 0; i < 10; i++){
        //Set random size and direction for each fish
        a[i] = (rand() % MAX_SIZE) + 1;
    }
    
    for(uint32_t i = 0; i < 10; i++){
        printf("%d,", a[i]);
    }
    puts("");

    int num = 2;
    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++){
            int less = -1;

            if(a[i] > a[j]) less = j;
            if(a[i] < a[j]) less = i;

            if(less >= 0){
                a[less] = a[num-1];
                j--;
                num--;
            }
        }
    }
    
    for(uint32_t i = 0; i < num; i++){
        printf("%d,", a[i]);
    }
    puts("");

    return 0;
}
