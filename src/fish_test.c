#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 100            //Number of fish (parameter)
#define L 10            //Dimension of fishtank (parameter)
#define v 1             //Speed of fish (parameter)
#define d 1             //Max eating distance (parameter)
#define t 2             //Time of each step (parameter)
#define MAX_SIZE 100
#define SECONDS_IN_DAY 2//86400 debug

typedef struct fish_s{
    int32_t x, y, z;
    uint16_t speed;
    uint16_t size;
    uint8_t direction;
}fish;

//Segment containing the local fish
fish *fish_local = NULL;
uint32_t num_fish_local = 0, length_fish_local = 0;
//Queue of fish that are passed to next local segment
fish *buffer_fish_out = NULL;
uint32_t index_out = 0, length_out = 0;
//Queue of fish that are passed by the previous segment
fish *buffer_fish_in = NULL;
uint32_t length_in = 0;
//Indexes of fish that are interacting for info on eaten fish
uint32_t *indexes_interaction = NULL;
uint32_t index_indexes = 0, length_indexes = 0;
//Max dimention of the z coordinate for the local segment (used as boudary between segments)
uint32_t max_z_local = 0;

fish* add_fish(fish *queue, uint32_t *index, uint32_t *length, fish toAdd){
    //If queue is full, increase size
    if((*index) >= (*length)){
        queue = (fish*)realloc(queue, (*length) * 2 * sizeof(fish));
        (*length) = (*length)*2;
    }

    queue[(*index)] = toAdd;
    (*index)++;

    return queue;
}

fish remove_fish(fish *queue, uint32_t *index, uint32_t toRemove){
    fish removed = queue[toRemove];
    queue[toRemove] = queue[(*index) - 1];
    (*index)--;

    return removed;
}

uint32_t* add_index(uint32_t *queue, uint32_t *index, uint32_t *length, uint32_t i){
    if((*index) >= (*length)){
        queue = (uint32_t*)realloc(queue, (*length) * 2 * sizeof(uint32_t));
        (*length) = (*length)*2;
    }

    queue[*index] = i;
    (*index)++;
    
    return queue;
}

double distance(fish f1, fish f2){
    int32_t x = f1.x - f2.x;
    int32_t y = f1.y - f2.y;
    int32_t z = f1.z - f2.z;

    return sqrt(x*x + y*y + z*z);
}

void eat_fish(fish *f1, fish *f2, uint32_t *i, uint32_t *j){
    fish *bigger = NULL;
    uint32_t smaller;
    
    if(f1->size > f2->size){
        bigger = f1;
        smaller = (*j);
    }
    else if(f1->size < f2->size){
        bigger = f2;
        smaller = (*i);
    }

    //If not same size the bigger increase size and the smaller is removed (he ded)
    if(bigger){
        bigger->size++;
        remove_fish(fish_local, &num_fish_local, smaller);
        (*j)--;                                                        //Necessary to confront all pairs (based on the way the remove works)
    }
}

void move_fish(uint32_t toMove){
    uint32_t distance = fish_local[toMove].speed * t;

    switch (fish_local[toMove].direction){
    
    case 0:
        fish_local[toMove].x = (fish_local[toMove].x + distance) % L;
        break;
    
    case 1:
        fish_local[toMove].y = (fish_local[toMove].y + distance) % L;
        break;
    
    default:
        //Fish would exit from local segment is put in the array of fish to send to next process
        if(fish_local[toMove].z + distance >= max_z_local){
            fish exiting_fish = remove_fish(fish_local, &num_fish_local, toMove);
            
            //Calculating new z coordinate in the next segment and add to exiting segment_local
            exiting_fish.z = (exiting_fish.z + distance) % max_z_local;
            buffer_fish_out = add_fish(buffer_fish_out, &index_out, &length_out, exiting_fish);
            break;
        }

        fish_local[toMove].z += distance;
        break;
    }
}






int main (int argc, char const** argv){
    MPI_Init(NULL, NULL);

    //Creating fish MPI datatype
    fish fish_definition;
    MPI_Datatype mpi_fish;
    int block_lens[6] = {1, 1, 1, 1, 1, 1};
    MPI_Datatype types[6] = {  MPI_INT32_T, MPI_INT32_T, MPI_INT32_T,
                                        MPI_UINT16_T, MPI_UINT16_T, MPI_UINT8_T};
    //Displacement
    MPI_Aint displacements[6];
    displacements[0] = offsetof(fish, x);
    displacements[1] = offsetof(fish, y);
    displacements[2] = offsetof(fish, z);
    displacements[3] = offsetof(fish, speed);
    displacements[4] = offsetof(fish, size);
    displacements[5] = offsetof(fish, direction);
    //Create mpi struct
    MPI_Type_create_struct(6, block_lens, displacements, types, &mpi_fish);
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
    
    //Allocation of local fish
    fish_local = (fish*)malloc(sizeof(fish)*num_fish_local);
    length_fish_local = num_fish_local;
    
    //Spawn of all the fishes of the local segment
    for(uint32_t i = 0; i < num_fish_local; i++){
        //Generate random coordinates and check if a fish is already present
        fish_local[i].x = rand() % L;
        fish_local[i].y = rand() % L;
        fish_local[i].z = rand() % max_z_local;

        //Set random size and direction for each fish
        fish_local[i].size = (rand() % MAX_SIZE) + 1;
        fish_local[i].direction = rand() % 3;
        fish_local[i].speed = v;
    }

    //Move the fish (also puts fish in fish_exiting)
    buffer_fish_out = (fish*)malloc(sizeof(fish) * (num_fish_local / 2));
    index_out = 0, length_out = num_fish_local;

    for (uint32_t i = 0; i < num_fish_local; i++){
        move_fish(i);
    }
    time_in_seconds += t;

    //Send the fish that are exiting the segment
    #define EXIT_FISH 32
    MPI_Request req1;
    MPI_Isend(buffer_fish_out, index_out, mpi_fish, (rank+1)%num_segments, EXIT_FISH, MPI_COMM_WORLD, &req1);
    //Receiving fish
    MPI_Status status_entering;
    MPI_Probe(MPI_ANY_SOURCE, EXIT_FISH, MPI_COMM_WORLD, &status_entering);
    MPI_Get_count(&status_entering, mpi_fish, &length_in);
    buffer_fish_in = malloc(sizeof(fish)*length_in);
    MPI_Recv(buffer_fish_in, length_in, mpi_fish, MPI_ANY_SOURCE, EXIT_FISH, MPI_COMM_WORLD, &status_entering);

    for(uint32_t i = 0; i < length_in; i++){
        fish_local = add_fish(fish_local, &num_fish_local, &length_fish_local, buffer_fish_in[i]);
    }
    free(buffer_fish_in);

    //Calculate the fish eaten locally
    if(num_fish_local >= 2){
        for(uint32_t i = 0; i < num_fish_local; i++){
            for(uint32_t j = i+1; j < num_fish_local; j++){
                fish *f1 = &fish_local[i];
                fish *f2 = &fish_local[j];
                if((int)distance(*f1, *f2) <= d){
                    eat_fish(f1, f2, &i, &j);
                }
            }
        }
    }
    
    //Checkign if the fish can interact with other in the next segment (distance from border < eating distance)
    index_out = 0;                                                      //Check the send is completed for safety TODO
    indexes_interaction = malloc(sizeof(fish) * length_out);
    index_indexes = 0, length_indexes = length_out;

    for(uint32_t i = 0; i < num_fish_local; i++){
        if(fish_local[i].z > max_z_local - d){
            //Adjusting coordinate z to allow next segment to calculate the distance correctly
            fish temp = fish_local[i];
            temp.z = temp.z - max_z_local;
            indexes_interaction = add_index(indexes_interaction, &index_indexes, &length_indexes, i);
            buffer_fish_out = add_fish(buffer_fish_out, &index_out, &length_out, temp);
        }
    }

    //Sending fish that can interact
    #define FISH_INTERACTION 33
    MPI_Request req2;
    MPI_Isend(buffer_fish_out, index_out, mpi_fish, (rank+1)%num_segments, FISH_INTERACTION, MPI_COMM_WORLD, &req2);
    //Receving fish that can interact
    MPI_Status status_interaction;
    MPI_Probe(MPI_ANY_SOURCE, FISH_INTERACTION, MPI_COMM_WORLD, &status_interaction);
    MPI_Get_count(&status_interaction, mpi_fish, &length_in);
    buffer_fish_in = malloc(sizeof(fish) * length_in);
    MPI_Recv(buffer_fish_in, length_in, mpi_fish, MPI_ANY_SOURCE, FISH_INTERACTION, MPI_COMM_WORLD, &status_interaction);

    //Calculate fish eaten in other segment
    #define EATEN 23
    for(uint32_t i = 0; i < length_in; i++){
        for(uint32_t j = 0; j < num_fish_local; j++){
            fish *f1 = &buffer_fish_in[i];
            fish *f2 = &fish_local[j];

            if(distance(*f1, *f2) <= d){
                //If f1 eats then remove f2 and return array of indexes of fish that have eaten
                if(f1->size > f2->size){
                    f1->size++;
                    remove_fish(fish_local, &num_fish_local, j);
                }
                //if f2 eats then set index on eaten fish array
                if(f1->size < f2->size){
                    f2->size++;
                    f1->direction = EATEN;
                }
            }
        }
    }
    //Preparing the buffer of fish to send back
    index_out = 0;                                  //Check previous send is completed for safety TODO
    for(uint32_t i = 0; i < length_in; i++){
        add_fish(buffer_fish_out, &index_out, &length_out, buffer_fish_in[i]);
    }
    free(buffer_fish_in);

    //Sending info of fish eaten to the segments
    if(rank != 0){
        MPI_Request req3;
        MPI_Isend(buffer_fish_out, index_out, mpi_fish, rank-1, FISH_INTERACTION, MPI_COMM_WORLD, &req3);      //Don't let last segment send stuff
    }
    if(rank != num_segments - 1){
        buffer_fish_in = malloc(sizeof(fish) * length_indexes);
        length_in = length_indexes;
        MPI_Recv(buffer_fish_in, length_in, mpi_fish, MPI_ANY_SOURCE, FISH_INTERACTION, MPI_COMM_WORLD, &status_interaction);
    }

    //Updating size and fish eaten according the the interaction result
    if(rank != num_segments - 1){
        for(uint32_t i = 0; i < index_indexes; i++){
            if(buffer_fish_in[i].direction = EATEN){
                remove_fish(fish_local, &num_fish_local, indexes_interaction[i]);
            }
            else{
                fish_local[indexes_interaction[i]].size = buffer_fish_in[i].size;
            }
        }
    }

    //If end of the day gather info on max and min size
    if(time_in_seconds >= SECONDS_IN_DAY){
        //Update day
        time_in_days++;
        time_in_seconds -= SECONDS_IN_DAY;

        //Calculate max and min size of the lake
        uint32_t max_size_local = 0;
        uint32_t min_size_local = 0xFF;
        //printf("Rank %d Num fish local %d, Segment length %d\n", rank, num_fish_local, length_segment_local);
        for(uint32_t i = 0; i < num_fish_local; i++){
            if(fish_local[i].size > max_size_local) max_size_local = fish_local[i].size;
            if(fish_local[i].size < min_size_local) min_size_local = fish_local[i].size;
            //printf("FIsh size: %d\n", segment_local[i].size);
        }
        uint16_t max_size_global;
        uint16_t min_size_global;
        MPI_Reduce(&max_size_local, &max_size_global, 1, MPI_UINT16_T, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&min_size_local, &min_size_global, 1, MPI_UINT16_T, MPI_MIN, 0, MPI_COMM_WORLD);

        if(rank == 0){
            printf("Day %d: Max size = %d, Min size = %d\n", time_in_days, max_size_global, min_size_global);
        }
    }
    
    //Free resources
    free(buffer_fish_out);
    free(fish_local);
    MPI_Type_free(&mpi_fish);

    MPI_Finalize();

    return 0;
}