CC=mpicc
FLAGS=-g -O3 -lm

all: hackathon

hackathon:
	$(CC) $(FLAGS) -o hackathon hackathon.c

clean:
	rm -f hackathon *~
