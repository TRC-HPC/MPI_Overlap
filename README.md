# MPI_Overlap
Build via:
make

Invoke using:
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1 -x UCX_TLS=self,sm,ud,rc ./overlap_regdata ...

usage:
./overlap_regdata -algorithm #data_elements #neighbors

where algorithm is one of:
-t: two-sided non-blocking comm. + wait for completion after compute (overlap),
-tw: two-sided comm. with immediate wait and then compute,
-o: one-sided asynchronous comm. + wait for completion after compute (overlap),
-ow: one-sided comm. with immediate wait and then compute

