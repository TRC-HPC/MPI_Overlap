# MPI_Overlap
Build via:
make

Invoke using:
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_TLS=self,sm,ud,rc ./overlap_regdata #commtype #data_elements #neighbors

where #commtype is one of:

  -t: two-sided non-blocking comm. + wait for completion after compute (overlap),

  -tw: two-sided comm. with immediate wait and then compute,

  -o: one-sided asynchronous comm. + wait for completion after compute (overlap),

  -ow: one-sided comm. with immediate wait and then compute

#data_elements is the number of data elements to send between ranks (typically ~200000-600000)

and #neighbors is the average number of neighbors each rank will send and receive data from.
