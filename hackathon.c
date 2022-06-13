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

#define rank0_printf(...) if  (nRank==0) {printf(__VA_ARGS__);}

// ---------------- //
// Global Variables //
// ---------------- //
// Number of elements in array
int g_nDataElementCount = 0;

// Number of send/receive neighbours
int g_nNeighbourCount = 0;

// World size
int g_nWorldSize = 0;
int g_nWorldHalfSize = 0;

// Communication Type
typedef enum
{
    COMM_ONESIDED_NONBLOCKING = 0,
    COMM_ONESIDED_BLOCKING,
    COMM_TWOSIDED_NONBLOCKING,
    COMM_TWOSIDED_BLOCKING,
    COMM_UNKNOWN,
} CommType_e;
CommType_e g_eCommType = COMM_UNKNOWN;


// Fill array with data
void fillArray(double* arr, int nValue)
{
    for (int i = 0; i < g_nDataElementCount; ++i)
	    arr[i] = nValue;
}

// Perform 'fake' computation to keep the cpu busy
void compute(double* array)
{
    // usleep(SLEEPTIME);
    // return;
    int nRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
    
    for (int  i = 0; i < COMP_ITRS; i++)
        for (int j = 0; j < g_nDataElementCount; j++)
            array[j] = nRank;
}

// Initialize one-sided window
void init_onesided(double* pExchangeArraySend, MPI_Win* window )
{
    MPI_Win_create(pExchangeArraySend, g_nDataElementCount*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, window);	
}

// Exchange one-sided data
void exchange_onesided(double* pExchangeArrayRecv,
        		       MPI_Request* req_arr, MPI_Win* window, int target)
{
    MPI_Rget(pExchangeArrayRecv, g_nDataElementCount, MPI_DOUBLE, target, 0, g_nDataElementCount, MPI_DOUBLE, *window, &req_arr[target]);
}

// Exchange two-sided nonblocking
void exchange_twosided(double* pExchangeArraySend, double* pExchangeArrayRecv, MPI_Request* req_arr, int nRank, int nTarget)
{
    const int nSendTag = MAX_MPI_NEIGHBORS * nRank + nTarget;		
    const int nRecvTag = nRank + MAX_MPI_NEIGHBORS * nTarget;

    MPI_Irecv(pExchangeArrayRecv, g_nDataElementCount, MPI_DOUBLE, nTarget, nRecvTag, MPI_COMM_WORLD, &req_arr[MAX_MPI_NEIGHBORS + nTarget]);
    MPI_Isend(pExchangeArraySend, g_nDataElementCount, MPI_DOUBLE, nTarget, nSendTag, MPI_COMM_WORLD, &req_arr[nTarget]);
}

// Extract commandline information
int parse_commandline(int argc, char** argv, int nRank)
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
        g_eCommType = COMM_ONESIDED_NONBLOCKING;
        rank0_printf("\tUsing Onesided communication\n");
    }
    else if (strcmp(argv[1], "-ow") == 0){
        g_eCommType = COMM_ONESIDED_BLOCKING;        
        rank0_printf("\tUsing Onesided_wait communication\n");
    }
    else if (strcmp(argv[1], "-t") == 0)
    {
        g_eCommType = COMM_TWOSIDED_NONBLOCKING;
        rank0_printf("\tUsing Twosided communication\n");
    }
    else if (strcmp(argv[1], "-tw") == 0)
    {
        g_eCommType = COMM_TWOSIDED_BLOCKING;
        rank0_printf("\tUsing Twosided_wait communication\n");
    }
    else
    {
        rank0_printf("Incorrect first argument given, choose from: onesided: -o, onesided_w: -ow, twosided: -t, twosided_w -tw\n");
        return -1;
    }
    
    g_nDataElementCount = atoi(argv[2]);
    if(g_nDataElementCount <= 0)
    {
        rank0_printf("Incorrect number of elements given: %d\n", g_nDataElementCount);
        return -1;
    }

    g_nNeighbourCount = atoi(argv[3]);
    if(g_nNeighbourCount < 2)
    {
        rank0_printf("Incorrect number of neighbors given: %d\n", g_nNeighbourCount);
        return -1;
    }

    return 0;
}


int do_communication_pattern(int *pTargetList, int nRank)
{
    // create target ranks pool for each rank
    int gtargets[g_nWorldSize][g_nWorldHalfSize];

    for(int k = 0; k < g_nWorldSize; k++)
        for(int j = 0; j < g_nWorldHalfSize; j++)
            gtargets[k][j] = -1;

    int nPossibleTargetList[g_nWorldHalfSize];

    // rank 0 determines the communication pattern
    // Each rank in the first group (0 .. g_nWorldHalfSize - 1) picks g_nNeighbourCount random ranks from the other group (g_nWorldHalfSize .. g_nWorldSize - 1)
    if(nRank == 0)
    {
        rank0_printf("Assigning communication:\n");
        srand(time(NULL));
        
        // Run on each rank of first group
        for(int crank = 0; crank < g_nWorldHalfSize; crank++)
        {
            // Initialize possible options
            for (int i = 0; i < g_nWorldHalfSize; i++)
                nPossibleTargetList[i] = i + g_nWorldHalfSize;

            // Try to pick g_nNeighbourCount different neighbours from second group
            for(int k = 0; k < g_nNeighbourCount; k++)
            {
                // Random Index
                const int nIndex = rand() % (g_nWorldHalfSize - k);

                // Get Target
                const int nTarget = nPossibleTargetList[nIndex];

                // Place the last index instead
                nPossibleTargetList[nIndex] = nPossibleTargetList[g_nWorldHalfSize - k - 1];
                    
                // Get Random Target
                gtargets[crank][k] = nTarget;
            }
        }
        

        // assign reverse targets to other group (world_half .. world_size - 1) based on symmetric communication pattern
        // Note: certain ranks in the second group will have to communicate with more than g_nNeighbourCount neighbors, due to random assignments of first group
        int nMaxNeighbourCount = 0;
        for(int crank = g_nWorldHalfSize; crank < g_nWorldSize; crank++)
        {
            int nNeighbourCount = 0;
            for(int target = 0; target < g_nWorldHalfSize; target++)
            {    
                for(int k = 0; k < g_nNeighbourCount; k++)
                {    
                    if(gtargets[target][k] == crank)
                        gtargets[crank][nNeighbourCount++] = target;
                }
            }

            // Update maximum neighbour count
            if(nNeighbourCount > nMaxNeighbourCount)
                nMaxNeighbourCount = nNeighbourCount;
        }
        rank0_printf("\tMaximum number of neighbors: %d\n", nMaxNeighbourCount);
    }

    rank0_printf("Distributing communication targets to worker ranks.\n");

    // distribute targets to individual ranks
    MPI_Scatter(gtargets , g_nWorldHalfSize , MPI_INT , pTargetList , g_nWorldHalfSize , MPI_INT , 0, MPI_COMM_WORLD);    

    return 0;
}

int main (int argc, char** argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get Ranks Info
    int nRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_nWorldSize);

    // Invalid World Size
    if (g_nWorldSize % 2)
    {
    	rank0_printf("Error: this program is meant to be run on an even number of ranks\n");
        MPI_Finalize();
        return 0;
    }

    // Half world size
    g_nWorldHalfSize = g_nWorldSize / 2;

    // Extract data from arguments
    if(parse_commandline(argc, argv, nRank))
    {
        MPI_Finalize();
        return 0;
    }

    // World half-size is not big enough to contain desired number of neighbours to communicate with    
    if (g_nWorldHalfSize < g_nNeighbourCount)
    {
    	rank0_printf("Error: g_nNeighbourCount=%d is too big compared to ranks/2 = %d\n", g_nNeighbourCount, g_nWorldHalfSize);
        MPI_Finalize();
        return 0;
    }

    // Initialize our send/recv targets
    int nTargetList[g_nWorldHalfSize];    
    if(do_communication_pattern(nTargetList, nRank))
    {
        MPI_Finalize();
        return 0;
    }

    // Init arrays
    rank0_printf("Allocating arrays:\n");
    rank0_printf("\tEach array (send/recv) is of size = %d bytes\n", sizeof(double)*g_nDataElementCount);

    // Compute array
    double *pComputeArray = (double *)malloc(g_nDataElementCount*sizeof(double));
    if(!pComputeArray)
    {
        printf("Error: unable to alloc arr");
        MPI_Finalize();
        return -1;
    }
    fillArray(pComputeArray, nRank);

    // Exchange Sent Array
    double *pExchangeArraySend = (double *)malloc(g_nDataElementCount*sizeof(double));
    if(!pExchangeArraySend)
    {
        printf("Error: unable to alloc pExchangeArraySend");
        MPI_Finalize();
        return -1;
    }
    fillArray(pExchangeArraySend, nRank);
    
    // Exchange Recv Arrays
    double *pExchangeArrayRecvList[g_nWorldSize];
    for(int k = 0; k < g_nWorldSize; k++)
    {
        pExchangeArrayRecvList[k] = (double *)malloc(g_nDataElementCount*sizeof(double));
        if(!pExchangeArrayRecvList[k])
        {
            printf("Error: unable to alloc pExchangeArrayRecvList[%u]", k);
            MPI_Finalize();
            return -1;
        }
        fillArray(pExchangeArrayRecvList[k], -1);
    }

    // MPI Requests
    MPI_Request* req_arr;
    int nreq = 2 * MAX_MPI_NEIGHBORS;
    req_arr = (MPI_Request*) malloc(nreq * sizeof(MPI_Request));
    if(!req_arr)
    {
        printf("Error: unable to alloc req_arr");
        MPI_Finalize();
        return -1;
    }
    // Initialize requests
    for(unsigned int k = 0; k < nreq; k++)
        req_arr[k] = MPI_REQUEST_NULL;

    // Average Iteration time
    double dAverageIterTimes[ITERATIONS];

    // Initiailize one-sided
    const bool bIsOneSided = (g_eCommType == COMM_ONESIDED_BLOCKING || g_eCommType == COMM_ONESIDED_NONBLOCKING);
    MPI_Win window;
    if (bIsOneSided)
	    init_onesided(pExchangeArraySend, &window);
    
    // Working variables
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
        if (bIsOneSided)
        {
        	MPI_Win_lock_all(MPI_MODE_NOCHECK, window);
            for(int k = 0; k < g_nWorldHalfSize; k++)
            {
                const int nTarget = nTargetList[k];
                if(nTarget == -1)
                    break;
                exchange_onesided(pExchangeArrayRecvList[nTarget], req_arr, &window, nTarget);
            }
            if (g_eCommType == COMM_ONESIDED_BLOCKING)
            {
                MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
                MPI_Win_unlock_all(window);
            }
        }
        else
        {
            // Two-Sided communication
            for(int k = 0; k < g_nWorldHalfSize; k++)
            {
                const int nTarget = nTargetList[k];
                if(nTarget == -1)
                    break;
                exchange_twosided(pExchangeArraySend, pExchangeArrayRecvList[nTarget], req_arr, nRank, nTarget);
            }

            if(g_eCommType == COMM_TWOSIDED_BLOCKING)
                MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
        }
        end = MPI_Wtime();
        send_time = end - start;
        
        // Compute 
        start = MPI_Wtime();
        compute(pComputeArray);
        end = MPI_Wtime();
        compute_time = end - start;

        // Waitall
        start = MPI_Wtime();
        if (g_eCommType == COMM_ONESIDED_NONBLOCKING)
        {
            MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
            MPI_Win_unlock_all(window);
        }
        else if (g_eCommType == COMM_TWOSIDED_NONBLOCKING) {
            MPI_Waitall(nreq, req_arr, MPI_STATUSES_IGNORE);
        }
        end = MPI_Wtime();
        wait_time = end - start;

        // Statistics 
        iter_end = MPI_Wtime();
        iter_time = iter_end - iter_start;

        // Check received data
#if 1
        int j;
        for (int k = 0; k < g_nWorldHalfSize; k++)
        {
            const int nTarget = nTargetList[k];
            if(nTarget == -1)
                break;

            for (j = 0; j < g_nDataElementCount; j++)
            {
                if(pExchangeArrayRecvList[nTarget][j] != nTarget)
                {
                    printf("rank %d, target %d: inconsitency in recv buffer at position %d = %f (supposed to be %f)\n", nRank, nTarget, j, pExchangeArrayRecvList[nTarget][j], nTarget);
                    break;
                }
            }
            // if (j == g_nDataElementCount)
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
        rank0_printf("\tAvg MPI Exchange Time: %.6e  +  Avg Compute Time: %.6e  +  Avg MPI Wait time: %.6e  ~ Avg Iteration Time: %.6e\n", global_send_sum/g_nWorldSize, global_compute_sum/g_nWorldSize, global_wait_sum/g_nWorldSize, global_iter_sum/g_nWorldSize);
        rank0_printf("\tMax MPI Exchange Time: %.6e  +  Max Compute Time: %.6e  +  Max MPI Wait time: %.6e  ~ Max Iteration Time: %.6e\n", max_send, max_compute, max_wait, max_iter);
        
        dAverageIterTimes[i] = global_iter_sum / g_nWorldSize;
        if(max_iter > global_max_iter && i >= WARMUP)
            global_max_iter = max_iter;
        
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Average iteration time
    double global_sum = 0;
    for (int i = WARMUP; i < ITERATIONS; ++i){
	    global_sum += dAverageIterTimes[i];
    }
    rank0_printf("Avg iteration time of %d iterations: %.6e\n", ITERATIONS, global_sum / (ITERATIONS - WARMUP));  
    rank0_printf("Max iteration time of %d iterations: %.6e\n", ITERATIONS, global_max_iter);  

    global_end = MPI_Wtime();
    rank0_printf("Total time elapsed: %.6e\n", global_end - global_start);
    
    // Cleanup
    if (bIsOneSided)
	    MPI_Win_free(&window);
  
    free(req_arr);
    free(pComputeArray);
    free(pExchangeArraySend);
    for(int k = 0; k < g_nWorldSize; k++)
        free(pExchangeArrayRecvList[k]);

    MPI_Finalize();
}
