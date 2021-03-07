CC = gcc
CFLAGS = -g -Wall -Wno-unused-variable -std=c99
Lib = -lm
EX = BEMT01

all : $(EX)
$(EX) : BEMT.o main.o
	$(CC) $(CFLAGS) -o $@ $^ $(Lib)

BEMT.o : BEMT.c
	$(CC) $(CFLAGS) -o $@ -c $<

main.o : main.c
	$(CC) $(CFLAGS) -o $@ -c $<

all:
	rm -rf ../bin && mkdir ../bin && mv $(EX) ../bin

clean:
	rm -rf *.o rm -rf ../bin
