CC=mpicc
FLAGS=-g -O3 -lm

all: regdata customdata

regdata: 
	$(CC) $(FLAGS) -o overlap_regdata overlap_regdata.c

customdata:
	$(CC) $(FLAGS) -o overlap_customdata overlap_customdata.c

clean:
	rm -f overlap *~
