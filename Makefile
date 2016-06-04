program=proj2

SRC=proj2.c

CC=gcc

CFLAGS = -std=c99 -Wall -pedantic -g -lm

${program}:
	${CC} ${SRC} -o ${program} ${CFLAGS}

clean:
	rm -f *.o ${program}
