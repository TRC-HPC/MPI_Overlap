CC=mpicc
FLAGS=-g -O3 -lm

all: regdata customdata hackathon

regdata: 
	$(CC) $(FLAGS) -o overlap_regdata overlap_regdata.c

hackathon:
	$(CC) $(FLAGS) -o hackathon hackathon.c

customdata:
	$(CC) $(FLAGS) -o overlap_customdata overlap_customdata.c

clean:
	rm -f overlap*. hackathon *~
