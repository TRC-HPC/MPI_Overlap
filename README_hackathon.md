# MPI_Overlap
Build via:
make

usage:
./overlap_regdata #commtype #data_elements #neighbors

where commtype is one of:

  -t: two-sided non-blocking comm. + wait for completion after compute (overlap),

  -tw: two-sided comm. with immediate wait and then compute,

  -o: one-sided asynchronous comm. + wait for completion after compute (overlap),

  -ow: one-sided comm. with immediate wait and then compute


Hackathon assignment #1:

Invoke overlap_regdata using two-sided communication, with UCX and the UD/RC protocols, for a buffer size of 400 kb and at least 6 nominal neighbors, e.g.:
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_TLS=self,sm,ud,rc ./overlap_regdata -t 50000 6

Demonstrate consistent avg and max iteration times of less than 0.1 seconds for at least three consecutive invokations of overlap_regdata with two-sided communication

Hackathon assignment #2:

Invoke overlap_regdata using one-sided communication, with UCX and the RC protocol, for a buffer size of 800 kb and at least 6 nominal neighbors, e.g.:
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_TLS=self,sm,rc ./overlap_regdata -o 100000 6

Demonstrate consistent avg and max iteration times of less than 0.12 seconds for at least three consecutive invokations of overlap_regdata with one-sided communication
(Bonus) demonstrate computation-communication overlap by comparing runs with -ow (immediate wait) and -o (delayed wait)

Hackathon assignment #3:

Invoke overlap_regdata using one-sided communication, with UCX and the RC protocol, for a buffer size of 8 Mb and at least 6 nominal neighbors, e.g.:
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_TLS=self,sm,rc ./overlap_regdata -o 1000000 6

Resolve ucx errors (Transport retry count exceeded) and demonstrate consistent iteration times for this scenario
