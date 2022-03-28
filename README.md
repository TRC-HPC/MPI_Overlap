# MPI_Overlap
Build via:
make

Invoke using:

twosided:
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1 -x UCX_TLS=self,sm,ud,rc ./overlap -t
twosided_w:
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1 -x UCX_TLS=self,sm,ud,rc ./overlap -tw
persistent:
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1 -x UCX_TLS=self,sm,ud,rc ./overlap -p
onesided:
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1  -x UCX_TLS=self,sm,ud,rc ./overlap -o
onesided_w:
	mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_NET_DEVICES=mlx5_0:1,mlx5_2:1  -x UCX_TLS=self,sm,ud,rc ./overlap -ow
