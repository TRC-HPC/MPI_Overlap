#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

#define SIZE       92400     // sizeof(cfdcell_v0)/sizeof(double)*400
#define SLEEPTIME  50000
#define ITERATIONS 300
#define WARMUP 5
#define COMP_ITRS 250
#define MAX_MPI_RANKS 		1000
#define MAX_MPI_NEIGHBORS 	1000
#define NGHB_SIZE           10

#define rank0_printf(...) if  (rank==0) {printf(__VA_ARGS__);}

void build_array(double* arr, int rank)
{
    for (int i = 0; i < SIZE; ++i)
    {
	    arr[i] = rank;
    }
}

void compute(double* array, MPI_Request* req_arr)
{
    // usleep(SLEEPTIME);
    // return;
    int flag = 0;
    double val = array[0];
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    for (int  i= 0; i < COMP_ITRS; ++i)
    {
        //MPI_Iprobe(target, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
        //MPI_Testall(2, req_arr, &flag, MPI_STATUSES_IGNORE);
        for (int j = 0; j < SIZE; ++j)
        {
            array[j] = rank;
        }
    }
}

void init_onesided(double* exchange_arr_send, MPI_Win* window )
{
    MPI_Win_create(exchange_arr_send, SIZE*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, window);	
}

void exchange_onesided(double* exchange_arr_send, double* exchange_arr_recv,
		       MPI_Request* req_arr, MPI_Win* window, int rank, int stride, int target, int world_size)
{
    int ret;
    {
        MPI_Rget(exchange_arr_recv, SIZE, MPI_DOUBLE, target, 0, SIZE, MPI_DOUBLE, *window, &req_arr[target]);
    }
}

void init_persistent(double* exchange_arr_send, double* exchange_arr_recv, MPI_Request* req_arr, int rank, int stride, int target, int world_size)
{
    int ret;
    MPI_Send_init(exchange_arr_send, SIZE, MPI_DOUBLE, target, 1, MPI_COMM_WORLD, &req_arr[0]);
    MPI_Recv_init(exchange_arr_recv, SIZE, MPI_DOUBLE, target, 1, MPI_COMM_WORLD, &req_arr[1]);
}

void exchange_persistent(MPI_Request* req_arr)
{
    MPI_Startall(2,req_arr);
}

void exchange_twosided_sync(double* exchange_arr_send, double* exchange_arr_recv, MPI_Request* req_arr, int rank, int stride, int target, int world_size)
{
    int ret;
    if (rank < stride)
    {
        MPI_Send(exchange_arr_send, SIZE, MPI_DOUBLE, target, 1, MPI_COMM_WORLD);
        MPI_Recv(exchange_arr_recv, SIZE, MPI_DOUBLE, target, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Recv(exchange_arr_recv, SIZE, MPI_DOUBLE, target, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(exchange_arr_send, SIZE, MPI_DOUBLE, target, 1, MPI_COMM_WORLD);
    }
}

void chunk_send(double* buffer, int count, int chunk, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request)
{
	int idx = 0;

	for(idx = 0; idx < count; idx += chunk)
		MPI_Isend( &buffer[idx], fmin(chunk, count - idx), datatype, dest, tag++, comm, request++);
}

void chunk_recv(double* buffer, int count, int chunk, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request)
{
	int size = sizeof(datatype), idx = 0;
	for(idx = 0; idx < count; idx += chunk)
		MPI_Irecv( &buffer[idx], fmin(chunk, count - idx), datatype, source, tag++, comm, request++);
}

void exchange_twosided(double* exchange_arr_send, double* exchange_arr_recv, MPI_Request* req_arr, int rank, int stride, int target, int world_size)
{
    int ret;

    int send_tag = MAX_MPI_NEIGHBORS * rank + target;		
    int recv_tag = rank + MAX_MPI_NEIGHBORS * target;

    MPI_Irecv(exchange_arr_recv, SIZE, MPI_DOUBLE, target, recv_tag, MPI_COMM_WORLD, &req_arr[MAX_MPI_NEIGHBORS + target]);
    MPI_Isend(exchange_arr_send, SIZE, MPI_DOUBLE, target, send_tag, MPI_COMM_WORLD, &req_arr[target]);
}

int main (int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, world_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size % 2)
    {
    	rank0_printf("Error: this program is meant to be run on an even number of ranks\n");
    	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    	MPI_Finalize();	    
    }

    int stride = world_size / 2;

    if (stride <= NGHB_SIZE)
    {
    	rank0_printf("Error: NGHB_SIZE=%d is too big compared to ranks/2 = %d\n", NGHB_SIZE, stride);
    	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    	MPI_Finalize();	    
    }

    rank0_printf("Each send is of size = %d\n", sizeof(double)*SIZE);

    /* create target ranks pool for each rank */
    int gtargets[world_size][stride];
    int targets[stride];

    for(int k = 0; k < world_size; k++)
        for(int j = 0; j < stride; j++)
            gtargets[k][j] = -1;


    if(rank == 0)
    {
        srand(time(NULL));
        for(int crank = 0; crank < stride; crank++)
        {
            for(int k = 0; k < NGHB_SIZE; k++)
            {    
                bool success = false;
                do
                {
                    gtargets[crank][k] = rand() % (stride) + stride;

                    int i;
                    for(i = 0; i < k; i++)
                        if(gtargets[crank][k] == gtargets[crank][i])      // avoid previously assigned rank
                            break;
                    if(i == k)
                        success = true;
                } while (success == false);   

                // rank0_printf("crank %d: target %d\n", crank, gtargets[crank][k]);
            }
        }

        // assign reverse targets
        for(int crank = stride; crank < world_size; crank++)
        {
            int nghb = 0;
            for(int target = 0; target < stride; target++)
            {    
                for(int k = 0; k < NGHB_SIZE; k++)
                {    
                    if(gtargets[target][k] == crank)
                        gtargets[crank][nghb++] = target;
                }
            }
        }
    }

    // distribute targets to individual ranks
    MPI_Scatter( gtargets , stride , MPI_INT , targets , stride , MPI_INT , 0, MPI_COMM_WORLD);    

    /* Init arrays */
    double *arr = (double *)malloc(SIZE*sizeof(double));
    double *exchange_arr_send = (double *)malloc(SIZE*sizeof(double));
    
    double *exchange_arr_recv[world_size];
    for(int k = 0; k < world_size; k++)
    {
        exchange_arr_recv[k] = (double *)malloc(SIZE*sizeof(double));
        if(!exchange_arr_recv[k])
        {
            printf("Error: unable to alloc exchange_arr_recv[%u]", k);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            MPI_Finalize();	    
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
    	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    	MPI_Finalize();	    
    }

    for(unsigned int k = 0; k < nreq; k++)
        req_arr[k] = MPI_REQUEST_NULL;

    double iter_times[ITERATIONS];

    /* Get configuration */
    bool onesided = false;
    bool onesided_w = false;
    bool twosided = false;
    bool twosided_w = false;
    bool persistent = false;
    if (argc <= 1)
    {
        printf("No argument given, choose from: onesided: -o, onesided_w: -ow, twosided: -t, twosided_w -tw");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        MPI_Finalize();	    
    }
    if (strcmp(argv[1], "-o") == 0)
    {
        onesided = true;
        rank0_printf("Using Onesided communication\n");
    }
    else if (strcmp(argv[1], "-ow") == 0){
        onesided_w = true;
        rank0_printf("Using Onesided_wait communication\n");
    }
    else if (strcmp(argv[1], "-t") == 0)
    {
        twosided = true;
        rank0_printf("Using Twosided communication\n");
    }
    else if (strcmp(argv[1], "-tw") == 0)
    {
        twosided_w = true;
        rank0_printf("Using Twosided_wait communication\n");
    }
    else if (strcmp(argv[1], "-p") == 0)
    {
        twosided_w = true;
        rank0_printf("Using Persistent communication\n");
    }    
    else
    {
        printf("Incorrect argument given, choose from: onesided: -o, onesided_w: -ow, twosided: -t, twosided_w -tw, persistent -p\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        MPI_Finalize();	    
    }
    
    /* Init onesided or persistent*/
    MPI_Win window;
    if (onesided || onesided_w)
    {
	    init_onesided(exchange_arr_send, &window);
    }
    else if (persistent)
    {
        for(int k = 0; k < NGHB_SIZE; k++)
        {
            int target = targets[k];
	        init_persistent(exchange_arr_send, exchange_arr_recv[target], req_arr, rank, stride, target, world_size);
        }
    }
    
    double start, end, iter_start, iter_end, compute_time, wait_time, send_time, iter_time;
    double global_send_sum, global_compute_sum, global_wait_sum, global_iter_sum;

    for (int i = 0; i < ITERATIONS; ++i)
    {
        iter_start = MPI_Wtime();
        if (rank == 0)
        {
            rank0_printf("Iteration %d\n", i);
        }
        
        /* Send/recv */
        start = MPI_Wtime();
        if (onesided || onesided_w)
        {
        	MPI_Win_lock_all(MPI_MODE_NOCHECK, window);
            for(int k = 0; k < stride; k++)
            {
                int target = targets[k];
                if(target == -1)
                    break;
                exchange_onesided(exchange_arr_send, exchange_arr_recv[target], req_arr, &window, rank, stride, target, world_size);
            }
            if (onesided_w)
            {
                MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
                MPI_Win_unlock_all(window);
            }
        }
        else if(twosided_w || twosided)
        {
            for(int k = 0; k < stride; k++)
            {
                int target = targets[k];
                if(target == -1)
                    break;
                exchange_twosided(exchange_arr_send, exchange_arr_recv[target], req_arr, rank, stride, target, world_size);
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
        compute(arr, req_arr);
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
        for (int k = 0; k < stride; k++)
        {
            int target = targets[k];
            if(target == -1)
                break;

            for (j = 0; j < SIZE; j++)
            {
                if(exchange_arr_recv[target][j] != target)
                {
                    printf("rank %d, target %d: inconsitency in recv buffer at position %d = %f (supposed to be %f)\n", rank, target, j, exchange_arr_recv[target][j], target);
                    break;
                }
            }
            // if (j == SIZE)
            //     printf("rank %d received all data correctly from target rank %d\n", rank, target);
        }
#endif

        MPI_Reduce(&send_time, &global_send_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&compute_time, &global_compute_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&wait_time, &global_wait_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&iter_time, &global_iter_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        rank0_printf("AVG - Exchange time: %f  |  Compute Time: %f  |  Wait time: %f  | Iteration Time: %f\n", global_send_sum/world_size, global_compute_sum/world_size, global_wait_sum/world_size, global_iter_sum/world_size);
        iter_times[i] = global_iter_sum/world_size;
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* Average iteration time */
    double global_sum = 0;
    for (int i = WARMUP; i < ITERATIONS; ++i){
	    global_sum += iter_times[i];
    }
    rank0_printf("Average iteration time: %f\n", global_sum / (ITERATIONS - WARMUP));  

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
