#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 20            //Number of fish (parameter)
#define L 10           //Dimension of fishtank (parameter)
#define v 1             //Speed of fish (parameter)
#define d 1             //Max eating distance (parameter)
#define t 2             //Time of each step (parameter)
#define MAX_SIZE 1000
#define SECONDS_IN_DAY 2//86400 debug

typedef struct fish{
    int32_t x, y, z;
    uint16_t speed;
    uint16_t size;
    uint8_t direction;
}fish;

//Segment containing the local fish
fish *segment_local = NULL;
uint32_t num_fish_local = 0, length_segment_local = 0;
//Queue of fish that are exiting the local segment
fish *fish_exiting = NULL;
uint32_t index_exit = 0, length_exit = 0;
//Queue of fish that are entering the local segment
fish *fish_entering = NULL;
uint32_t length_entering = 0;
//Queue of fish that can interact with fish from the next segment
fish *fish_interaction_out = NULL;
uint32_t index_interaction = 0, length_interaction = 0;
//Queue of fish from the previous segment that can interact with the local
fish *fish_interaction_in = NULL;
uint32_t length_interaction_in = 0;
//Max dimention of the z coordinate for the local segment (used as boudary between segments)
uint32_t max_z_local = 0;

void add_fish(fish *queue, uint32_t *index, uint32_t *length, fish toAdd){
    //If queue is full, increase size
    if(index == length){
        queue = (fish*)realloc(queue, (*length)*3/2);
    }

    queue[(*index)] = toAdd;
    (*index)++;
}

fish remove_fish(fish *queue, uint32_t *index, uint32_t toRemove){
    fish removed = queue[toRemove];
    queue[toRemove] = queue[(*index) - 1];
    (*index)--;

    return removed;
}

double distance(fish f1, fish f2){
    int32_t x = f1.x - f2.x;
    int32_t y = f1.y - f2.y;
    int32_t z = f1.z - f2.z;

    return sqrt(x*x + y*y + z*z);
}

void move_fish(uint32_t toMove){
    uint32_t distance = segment_local[toMove].speed * t;

    switch (segment_local[toMove].direction){
    
    case 0:
        segment_local[toMove].x = (segment_local[toMove].x + distance) % L;
        break;
    
    case 1:
        segment_local[toMove].y = (segment_local[toMove].y + distance) % L;
        break;
    
    default:
        //Fish would exit from local segment is put in the array of fish to send to next process
        if(segment_local[toMove].z + distance >= max_z_local){
            fish exiting_fish = remove_fish(segment_local, &num_fish_local, toMove);
            
            //Calculating new z coordinate in the next segment and add to exiting segment_local
            exiting_fish.z = (exiting_fish.z + distance) % max_z_local;
            add_fish(fish_exiting, &index_exit, &length_exit, exiting_fish);
            break;
        }

        segment_local[toMove].z += distance;
        
        //Fish can be eat by another fish in the near segment
        if(segment_local[toMove].z > max_z_local - d){
            //Adjusting coordinate z to allow next segment to calculate the distance correctly
            fish temp = segment_local[toMove];
            temp.z = temp.z - max_z_local;
            add_fish(fish_interaction_out, &index_interaction, &length_interaction, temp);
        }
        break;
    }
}

int main (int argc, char** argv){
    MPI_Init(NULL, NULL);

    //Creating fish MPI datatype
    fish fish_definition;
    MPI_Datatype mpi_fish;
    int struct_len = 6;
    int block_lens[struct_len];
    MPI_Datatype types[struct_len];
    //Displacement
    MPI_Aint displacements[struct_len];
    MPI_Aint current_displacement = 0;
    //3 coordinates
    block_lens[0] = 1;
    types[0] = MPI_INT32_T;
    displacements[0] = (size_t) &(fish_definition.x) - (size_t) &fish_definition;
    block_lens[1] = 1;
    types[1] = MPI_INT32_T;
    displacements[1] = (size_t) &(fish_definition.y) - (size_t) &fish_definition;
    block_lens[2] = 1;
    types[2] = MPI_INT32_T;
    displacements[2] = (size_t) &(fish_definition.z) - (size_t) &fish_definition;
    //Speed and size
    block_lens[3] = 1;
    types[3] = MPI_UINT16_T;
    displacements[3] = (size_t) &(fish_definition.speed) - (size_t) &fish_definition;
    block_lens[4] = 1;
    types[4] = MPI_UINT16_T;
    displacements[4] = (size_t) &(fish_definition.size) - (size_t) &fish_definition;
    //Direction
    block_lens[5] = 1;
    types[5] = MPI_UINT8_T;
    displacements[5] = (size_t) &(fish_definition.direction) - (size_t) &fish_definition;
    //Create mpi struct
    MPI_Type_create_struct(struct_len, block_lens, displacements, types, &mpi_fish);
    MPI_Type_commit(&mpi_fish);
    
    //Obtaining number fo processes and rank
    int rank;
    int num_segments;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_segments);
    srand(time(NULL) + rank);

    //Set local variables
    //Z dim for each segment in meters
    max_z_local = L / num_segments;
    if(rank == 0){
        max_z_local += L % num_segments;
    }

    //Num of fish for each segment (used as index of effective length of segment_local)
    num_fish_local = N / num_segments;
    if(rank == 0){
        num_fish_local += N % num_segments;
    }

    //Time
    uint32_t time_in_seconds = 0;
    uint32_t time_in_days = 0;
    
    //Allocation of queues
    segment_local = (fish*)malloc(sizeof(fish)*num_fish_local);
    length_segment_local = num_fish_local;

    fish_exiting = (fish*)malloc(sizeof(fish)*num_fish_local);
    index_exit = 0, length_exit = num_fish_local;

    fish_interaction_out = (fish*)malloc(sizeof(fish)*num_fish_local);
    index_interaction = 0, length_interaction = num_fish_local;
    
    //Initialization
    memset(segment_local, 0, sizeof(fish)*num_fish_local);
    memset(fish_exiting, 0, sizeof(fish)*num_fish_local);
    memset(fish_interaction_out, 0, sizeof(fish)*num_fish_local);
    
    uint16_t max_size_local = 0, min_size_local = 0xFF;

    //Spawn of all the fishes of the local segment
    for(uint32_t i = 0; i < num_fish_local; i++){
        //Generate random coordinates and check if a fish is already present
        segment_local[i].x = rand() % L;
        segment_local[i].y = rand() % L;
        segment_local[i].z = rand() % max_z_local;

        //Set random size and direction for each fish
        segment_local[i].size = (rand() % MAX_SIZE) + 1;
        segment_local[i].direction = rand() % 3;
        segment_local[i].speed = v;

        //Updating local min and max size
        if(segment_local[i].size > max_size_local){
            max_size_local = segment_local[i].size;
        }
        if(segment_local[i].size < min_size_local){
            min_size_local = segment_local[i].size;
        }
    }
    printf("Rank %d Num local fish: %d\n", rank, num_fish_local);

    //Move the fish
    for (uint32_t i = 0; i < num_fish_local; i++){
        move_fish(i);
    }
    time_in_seconds += t;

    //Send the fishes that are exiting the segment

    //Send the fish that are on the border of the segment (distance from border < eating distance)

    //Calculate the fish eaten locally (check there are at least 2 fish TODO)
    if(num_fish_local >= 2){
        for(uint32_t i = 0; i < num_fish_local; i++){
            for(uint32_t j = i+1; j < num_fish_local; j++){
                if((int)distance(segment_local[i], segment_local[j]) <= d){
                    fish *bigger = NULL;
                    uint32_t smaller;
                    if(segment_local[i].size > segment_local[j].size){
                        bigger = &segment_local[i];
                        smaller = j;
                    }
                    else if(segment_local[i].size < segment_local[j].size){
                        bigger = &segment_local[j];
                        smaller = i;
                    }

                    //If not same size the bigger increase size and the smaller is removed (he ded)
                    if(bigger){
                        bigger->size++;
                        remove_fish(segment_local, &num_fish_local, smaller);
                        j--;
                        printf("Rank %d, new size %d\n", rank, bigger->size);
                    }
                }
            }
        }
    }
    printf("Rank %d Num local fish: %d\n", rank, num_fish_local);

    //Calculate fish eaten in other segment

    //Sending info of fish eaten to the segments

    //If end of the day gather info on max and min size
    if(time_in_seconds >= SECONDS_IN_DAY){
        //Update day
        time_in_days++;
        time_in_seconds -= SECONDS_IN_DAY;

        //Calculate max and min size of the lake
        max_size_local = 0;
        min_size_local = 0xFF;
        for(uint32_t i = 0; i < num_fish_local; i++){
            if(segment_local[i].size > max_size_local) max_size_local = segment_local[i].size;
            if(segment_local[i].size < min_size_local) min_size_local = segment_local[i].size;
        }
        uint16_t max_size_global;
        uint16_t min_size_global;
        MPI_Reduce(&max_size_local, &max_size_global, 1, MPI_UINT16_T, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&min_size_local, &min_size_global, 1, MPI_UINT16_T, MPI_MIN, 0, MPI_COMM_WORLD);

        if(rank == 0){
            printf("Day %d: Max size = %d, Min size = %d\n", time_in_days, max_size_global, min_size_global);
        }
    }
    
    MPI_Type_free(&mpi_fish);
    MPI_Finalize();

    return 0;
}