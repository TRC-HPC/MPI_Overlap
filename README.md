# MPI_Overlap
Build via:
make

Invoke using:

twosided (wait after compute):
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1 -x UCX_TLS=self,sm,ud,rc ./overlap -t

twosided (immediate wait):
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1 -x UCX_TLS=self,sm,ud,rc ./overlap -tw

onesided (wait after compute):
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1  -x UCX_TLS=self,sm,ud,rc ./overlap -o

onesided (immediate wait):
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1  -x UCX_TLS=self,sm,ud,rc ./overlap -ow

