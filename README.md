# HPC Hackathon

## Prerequisites:
The following libraries are required:
  - [UCX](https://github.com/openucx/ucx)
  - [Open MPI](https://github.com/open-mpi/ompi.git) compiled with UCX support (./configure --prefix=*ompi-install-path* --with-ucx=*ucx-install-path*)

Alternatively, you can use the pre-built UCX and OpenMPI libraries by adding the following to your ~/.bashrc:
- Team 1:
```
HACK_DIR=/mnt/beegfs/hackathon/team1
export PATH=${HACK_DIR}/ucx/bin:${HACK_DIR}/ompi/bin:${PATH}
export LD_LIBRARY_PATH=${HACK_DIR}/ucx/lib:${HACK_DIR}/ompi/lib:${LD_LIBRARY_PATH}
```

- Team 2:
```
HACK_DIR=/mnt/beegfs/hackathon/team2
export PATH=${HACK_DIR}/ucx/bin:${HACK_DIR}/ompi/bin:${PATH}
export LD_LIBRARY_PATH=${HACK_DIR}/ucx/lib:${HACK_DIR}/ompi/lib:${LD_LIBRARY_PATH}
```

- Team 3:
```
HACK_DIR=/mnt/central/hackathon/bin
export PATH=${HACK_DIR}/ucx/bin:${HACK_DIR}/ompi/bin:${PATH}
export LD_LIBRARY_PATH=${HACK_DIR}/ucx/lib:${HACK_DIR}/ompi/lib:${LD_LIBRARY_PATH}
```

## Obtaining the code
git clone https://github.com/TRC-HPC/MPI_Overlap -b hackathon

## Running the code
Build via:
```
make
```

usage:
```
./hackathon #commtype #data_elements #neighbors
```

where commtype is one of:

  -t: two-sided non-blocking comm. + wait for completion after compute (overlap),

  -tw: two-sided comm. with immediate wait and then compute,

  -o: one-sided asynchronous comm. + wait for completion after compute (overlap),

  -ow: one-sided comm. with immediate wait and then compute

## Hackathon assignment #1:

Invoke hackathon on 5 nodes using two-sided communication, with UCX and the UD/RC protocols, for a buffer size of 400 kb and at least 6 nominal neighbors, e.g.:

```
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_TLS=self,sm,ud,rc ./hackathon -t 50000 6
```

Demonstrate consistent avg and max iteration times of less than 0.1 seconds for at least three consecutive invokations of hackathon with two-sided communication

## Hackathon assignment #2:

Invoke hackathon on 5 nodes using one-sided communication, with UCX and the RC protocol, for a buffer size of 800 kb and at least 6 nominal neighbors, e.g.:

```
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_TLS=self,sm,rc ./hackathon -o 100000 6
```

Demonstrate consistent avg and max iteration times of less than 0.12 seconds for at least three consecutive invokations of hackathon with one-sided communication

(Bonus) demonstrate computation-communication overlap by comparing runs with -ow (immediate wait) and -o (delayed wait)

## Hackathon assignment #3:

Invoke hackathon on 5 nodes using one-sided communication, with UCX and the RC protocol, for a buffer size of 8 Mb and at least 6 nominal neighbors, e.g.:

```
mpirun --hostfile hostfile --bind-to core --mca pml ucx -x UCX_IB_GID_INDEX=0 -x UCX_TLS=self,sm,rc ./hackathon -o 1000000 6
```

Resolve ucx errors (Transport retry count exceeded) and demonstrate consistent iteration times for this scenario

## Debugging (Advanced)
To debug the code, first make sure to change the FLAGS in Makefile to:
```
FLAGS="-g -O0 -lm" to disable optimizations.
```
Then, build debug versions of ucx and openmpi, or change the path in your ~/.bashrc to the pre-built development versions:
```
export PATH=${HACK_DIR}/ucx-devel/bin:${HACK_DIR}/ompi-debug/bin:${PATH}
export LD_LIBRARY_PATH=${HACK_DIR}/ucx-devel/lib:${HACK_DIR}/ompi-debug/lib:${LD_LIBRARY_PATH}
```
Debugging can be accomplished by attaching to one (or all) of the MPI processes using gdb, as described in:
(https://www.open-mpi.org/faq/?category=debugging)
Or by using VS code on a remote target machine, as shown in:
(https://iamsorush.com/posts/debug-mpi-vs-code/)
