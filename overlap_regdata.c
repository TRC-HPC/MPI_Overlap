#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

#define SLEEPTIME  50000		// usec
#define ITER_TRSH  1.0
#define ITERATIONS 200
#define WARMUP 5
#define COMP_ITRS 3000
#define MAX_MPI_RANKS 		1000
#define MAX_MPI_NEIGHBORS 	1000

#define rank0_printf(...) if  (rank==0) {printf(__VA_ARGS__);}

int data_elements = 0, nghb_size = 0;
int world_size = 0, world_half = 0;
bool onesided = false;
bool onesided_w = false;
bool twosided = false;
bool twosided_w = false;
bool persistent = false;


void build_array(double* arr, int rank)
{
    for (int i = 0; i < data_elements; ++i)
    {
	    arr[i] = rank;
    }
}

void compute(int* targets, double* array, MPI_Request* req_arr)
{
    // usleep(SLEEPTIME);
    // return;

    int flag = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    for (int  i = 0; i < COMP_ITRS; i++)
    {
        for (int j = 0; j < data_elements; j++)
        {
            array[j] = rank;
        }
    }
}

void init_onesided(double* exchange_arr_send, MPI_Win* window )
{
    MPI_Win_create(exchange_arr_send, data_elements*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, window);	
}

void exchange_onesided(double* exchange_arr_send, double* exchange_arr_recv,
        		       MPI_Request* req_arr, MPI_Win* window, int rank, int world_half, int target, int world_size)
{
    int ret;
    {
        MPI_Rget(exchange_arr_recv, data_elements, MPI_DOUBLE, target, 0, data_elements, MPI_DOUBLE, *window, &req_arr[target]);
    }
}

void init_persistent(double* exchange_arr_send, double* exchange_arr_recv, MPI_Request* req_arr, int rank, int world_half, int target, int world_size)
{
    int ret;
    MPI_Send_init(exchange_arr_send, data_elements, MPI_DOUBLE, target, 1, MPI_COMM_WORLD, &req_arr[0]);
    MPI_Recv_init(exchange_arr_recv, data_elements, MPI_DOUBLE, target, 1, MPI_COMM_WORLD, &req_arr[1]);
}

void exchange_persistent(MPI_Request* req_arr)
{
    MPI_Startall(2,req_arr);
}

void exchange_twosided_sync(double* exchange_arr_send, double* exchange_arr_recv, MPI_Request* req_arr, int rank, int world_half, int target, int world_size)
{
    int ret;
    int send_tag = MAX_MPI_NEIGHBORS * rank + target;		
    int recv_tag = rank + MAX_MPI_NEIGHBORS * target;

    if (rank < world_half)
    {
        MPI_Send(exchange_arr_send, data_elements, MPI_DOUBLE, target, send_tag, MPI_COMM_WORLD);
        MPI_Recv(exchange_arr_recv, data_elements, MPI_DOUBLE, target, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Recv(exchange_arr_recv, data_elements, MPI_DOUBLE, target, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(exchange_arr_send, data_elements, MPI_DOUBLE, target, send_tag, MPI_COMM_WORLD);
    }
}

void exchange_twosided(double* exchange_arr_send, double* exchange_arr_recv, MPI_Request* req_arr, int rank, int world_half, int target, int world_size)
{
    int ret;

    int send_tag = MAX_MPI_NEIGHBORS * rank + target;		
    int recv_tag = rank + MAX_MPI_NEIGHBORS * target;

    MPI_Irecv(exchange_arr_recv, data_elements, MPI_DOUBLE, target, recv_tag, MPI_COMM_WORLD, &req_arr[MAX_MPI_NEIGHBORS + target]);
    MPI_Isend(exchange_arr_send, data_elements, MPI_DOUBLE, target, send_tag, MPI_COMM_WORLD, &req_arr[target]);
}

int parse_commandline(int argc, char** argv, int rank)
{
    /* parse commandline arguments */
    rank0_printf("Parsing commandline arguments:\n");
    if (argc < 4)
    {
        rank0_printf("Not enough arguments given.\n");
        rank0_printf("Usage: %s %s %s %s\n", argv[0], "-algorithm (o,ow,t,tw)", "#data_elements", "#neighbors");
        return -1;
    }
    if (strcmp(argv[1], "-o") == 0)
    {
        onesided = true;
        rank0_printf("\tUsing Onesided communication\n");
    }
    else if (strcmp(argv[1], "-ow") == 0){
        onesided_w = true;
        rank0_printf("\tUsing Onesided_wait communication\n");
    }
    else if (strcmp(argv[1], "-t") == 0)
    {
        twosided = true;
        rank0_printf("\tUsing Twosided communication\n");
    }
    else if (strcmp(argv[1], "-tw") == 0)
    {
        twosided_w = true;
        rank0_printf("\tUsing Twosided_wait communication\n");
    }
    else if (strcmp(argv[1], "-p") == 0)
    {
        twosided_w = true;
        rank0_printf("\tUsing Persistent communication\n");
    }    
    else
    {
        rank0_printf("Incorrect first argument given, choose from: onesided: -o, onesided_w: -ow, twosided: -t, twosided_w -tw, persistent -p\n");
        return -1;
    }
    
    data_elements = atoi(argv[2]);
    if(data_elements <= 0)
    {
        rank0_printf("Incorrect number of elements given: %d\n", data_elements);
        return -1;
    }

    nghb_size = atoi(argv[3]);
    if(nghb_size < 2)
    {
        rank0_printf("Incorrect number of neighbors given: %d\n", nghb_size);
        return -1;
    }

    return 0;
}

int do_communication_pattern(int *targets, int world_size, int world_half, int rank)
{
    /* create target ranks pool for each rank */
    int gtargets[world_size][world_half];

    for(int k = 0; k < world_size; k++)
        for(int j = 0; j < world_half; j++)
            gtargets[k][j] = -1;


    /* rank 0 determines the communication pattern */
    if(rank == 0)
    {
        rank0_printf("Assigning communication targets:\n");
        srand(time(NULL));
        /* each rank in the first group (0 .. world_half - 1) picks nghb_size random ranks from the other group (world_half .. world_size - 1) */
        for(int crank = 0; crank < world_half; crank++)
        {
            for(int k = 0; k < nghb_size; k++)
            {    
                bool success = false;
                do
                {
                    gtargets[crank][k] = rand() % (world_half) + world_half;

                    int i;
                    for(i = 0; i < k; i++)
                        if(gtargets[crank][k] == gtargets[crank][i])      // avoid previously assigned rank
                            break;
                    if(i == k)
                        success = true;
                } while (success == false);   
            }
        }

        /* assign reverse targets to other group (world_half .. world_size - 1) based on symmetric communication pattern */
        /* note: certain ranks will have to communicate with more than nghb_size neighbors, due to random assignments of first group */
        int max_nghb = 0;
        for(int crank = world_half; crank < world_size; crank++)
        {
            int nghb = 0;
            for(int target = 0; target < world_half; target++)
            {    
                for(int k = 0; k < nghb_size; k++)
                {    
                    if(gtargets[target][k] == crank)
                        gtargets[crank][nghb++] = target;
                }
            }
            if(nghb > max_nghb)
                max_nghb = nghb;
        }
        rank0_printf("\tMaximum number of neighbors: %d\n", max_nghb);
    }

    rank0_printf("Distributing communication targets to worker ranks.\n");
    // distribute targets to individual ranks
    MPI_Scatter( gtargets , world_half , MPI_INT , targets , world_half , MPI_INT , 0, MPI_COMM_WORLD);    


    // output targets to file
#if 0
    char ranks_fname[50];
    sprintf(ranks_fname, "rank%d_targets.txt", rank);
    FILE *franks = fopen(ranks_fname, "w");
    if(franks == NULL)
    {
        printf("Error: unable to open file %s", ranks_fname);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        MPI_Finalize();	    
    }
    for(int k = 0; k < world_half; k++)
    {
        int target = targets[k];
        if(target == -1)
            break;
        fprintf(franks, "%d ", target);
    }
    fclose(franks);
#endif

    return 0;
}

int main (int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, world_size, world_half;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size % 2)
    {
    	rank0_printf("Error: this program is meant to be run on an even number of ranks\n");
        MPI_Finalize();
        return 0;
    }

    world_half = world_size / 2;

    if(parse_commandline(argc, argv, rank))
    {
        MPI_Finalize();
        return 0;
    }

    if (world_half <= nghb_size)
    {
    	rank0_printf("Error: nghb_size=%d is too big compared to ranks/2 = %d\n", nghb_size, world_half);
        MPI_Finalize();
        return 0;
    }

    int targets[world_half];
    
    if(do_communication_pattern(targets, world_size, world_half, rank))
    {
        MPI_Finalize();
        return 0;
    }

    /* Init arrays */
    rank0_printf("Allocating arrays:\n");
    rank0_printf("\tEach array (send/recv) is of size = %d bytes\n", sizeof(double)*data_elements);

    double *arr = (double *)malloc(data_elements*sizeof(double));
    if(!arr)
    {
        printf("Error: unable to alloc arr");
        MPI_Finalize();
        return -1;
    }

    double *exchange_arr_send = (double *)malloc(data_elements*sizeof(double));
    if(!exchange_arr_send)
    {
        printf("Error: unable to alloc exchange_arr_send");
        MPI_Finalize();
        return -1;
    }
    
    double *exchange_arr_recv[world_size];
    for(int k = 0; k < world_size; k++)
    {
        exchange_arr_recv[k] = (double *)malloc(data_elements*sizeof(double));
        if(!exchange_arr_recv[k])
        {
            printf("Error: unable to alloc exchange_arr_recv[%u]", k);
            MPI_Finalize();
            return -1;
        }
        build_array(exchange_arr_recv[k], -1);
    }

    build_array(arr, rank);
    build_array(exchange_arr_send, rank);

    MPI_Request* req_arr;
    int nreq = 2 * MAX_MPI_NEIGHBORS;
    req_arr = (MPI_Request*) malloc(nreq * sizeof(MPI_Request));
    if(!req_arr)
    {
        printf("Error: unable to alloc req_arr");
        MPI_Finalize();
        return -1;
    }

    for(unsigned int k = 0; k < nreq; k++)
        req_arr[k] = MPI_REQUEST_NULL;


    double avg_iter_times[ITERATIONS];

    /* Init onesided or persistent*/
    MPI_Win window;
    if (onesided || onesided_w)
    {
	    init_onesided(exchange_arr_send, &window);
    }
    else if (persistent)
    {
        for(int k = 0; k < nghb_size; k++)
        {
            int target = targets[k];
	        init_persistent(exchange_arr_send, exchange_arr_recv[target], req_arr, rank, world_half, target, world_size);
        }
    }
    
    double start, end, iter_start, iter_end, compute_time, wait_time, send_time, iter_time, global_start, global_end;
    double global_send_sum, global_compute_sum, global_wait_sum, global_iter_sum, global_max_iter = 0;
    double max_send, max_compute, max_wait, max_iter;

    rank0_printf("\n\nBeginning send/recv:\n");
    global_start = MPI_Wtime();

    for (int i = 0; i < ITERATIONS; ++i)
    {
        iter_start = MPI_Wtime();
        rank0_printf("Iteration %d\n", i);
        
        /* Send/recv */
        start = MPI_Wtime();
        if (onesided || onesided_w)
        {
        	MPI_Win_lock_all(MPI_MODE_NOCHECK, window);
            for(int k = 0; k < world_half; k++)
            {
                int target = targets[k];
                if(target == -1)
                    break;
                exchange_onesided(exchange_arr_send, exchange_arr_recv[target], req_arr, &window, rank, world_half, target, world_size);
            }
            if (onesided_w)
            {
                MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
                MPI_Win_unlock_all(window);
            }
        }
        else if(twosided_w || twosided)
        {
            for(int k = 0; k < world_half; k++)
            {
                int target = targets[k];
                if(target == -1)
                    break;
                exchange_twosided(exchange_arr_send, exchange_arr_recv[target], req_arr, rank, world_half, target, world_size);
            }

            if(twosided_w)
                MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
        }
        else if (persistent)
        {
            exchange_persistent(req_arr);
        }
        end = MPI_Wtime();
        send_time = end - start;
        
        /* Compute */
        start = MPI_Wtime();
        compute(targets, arr, req_arr);
        end = MPI_Wtime();
        compute_time = end - start;

        /* Waitall */
        start = MPI_Wtime();
        if (onesided)
        {
            MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
            MPI_Win_unlock_all(window);
        }
        else if (twosided || persistent) {
            MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
        }
        end = MPI_Wtime();
        wait_time = end - start;

        /* Statistics */
        iter_end = MPI_Wtime();
        iter_time = iter_end - iter_start;

        /* Check received data */
#if 1
        int j;
        for (int k = 0; k < world_half; k++)
        {
            int target = targets[k];
            if(target == -1)
                break;

            for (j = 0; j < data_elements; j++)
            {
                if(exchange_arr_recv[target][j] != target)
                {
                    printf("rank %d, target %d: inconsitency in recv buffer at position %d = %f (supposed to be %f)\n", rank, target, j, exchange_arr_recv[target][j], target);
                    break;
                }
            }
            // if (j == data_elements)
            //     printf("rank %d received all data correctly from target rank %d\n", rank, target);
        }
#endif

        MPI_Reduce(&send_time, &global_send_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&compute_time, &global_compute_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&wait_time, &global_wait_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&iter_time, &global_iter_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&send_time, &max_send, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&compute_time, &max_compute, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&wait_time, &max_wait, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&iter_time, &max_iter, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        rank0_printf("\tAvg MPI Exchange Time: %.6e  +  Avg Compute Time: %.6e  +  Avg MPI Wait time: %.6e  ~ Avg Iteration Time: %.6e\n", global_send_sum/world_size, global_compute_sum/world_size, global_wait_sum/world_size, global_iter_sum/world_size);
        rank0_printf("\tMax MPI Exchange Time: %.6e  +  Max Compute Time: %.6e  +  Max MPI Wait time: %.6e  ~ Max Iteration Time: %.6e\n", max_send, max_compute, max_wait, max_iter);
        
        avg_iter_times[i] = global_iter_sum / world_size;
        if(max_iter > global_max_iter && i >= WARMUP)
            global_max_iter = max_iter;
        
#if 0
        if(iter_time > ITER_TRSH)
             printf("rank %d: iter_time %f > %f\n", rank, iter_time, ITER_TRSH);
#endif

        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* Average iteration time */
    double global_sum = 0;
    for (int i = WARMUP; i < ITERATIONS; ++i){
	    global_sum += avg_iter_times[i];
    }
    rank0_printf("Avg iteration time of %d iterations: %.6e\n", ITERATIONS, global_sum / (ITERATIONS - WARMUP));  
    rank0_printf("Max iteration time of %d iterations: %.6e\n", ITERATIONS, global_max_iter);  

    global_end = MPI_Wtime();
    rank0_printf("Total time elapsed: %.6e\n", global_end - global_start);
    
    /* Cleanup */
    if (onesided || onesided_w){
	    MPI_Win_free(&window);
    }
  
    free(req_arr);
    free(arr);
    free(exchange_arr_send);
    for(int k = 0; k < world_size; k++)
        free(exchange_arr_recv[k]);

    MPI_Finalize();
}
