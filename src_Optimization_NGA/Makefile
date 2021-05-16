##
# Makefile for the FGA stuff
##
CC	= gcc
CFLAGS	= -O3
TARGETS	= Optimization_NGA
LIBS	= -lm
OBJS	= random.o forward.o ga_oper.o

Optimization_NGA:	fga.c ga_oper.h $(OBJS)
	$(CC) $(CFLAGS) -o $@ fga.c $(OBJS) $(LIBS)

ga_oper.o: ga_oper.h random.o forward.o
	$(CC) $(CFLAGS) -c ga_oper.c

random.o: random.h
	$(CC) $(CFLAGS) -c random.c

forward.o: forward.h $(FOROBJ)
	$(CC) $(CFLAGS) -c forward.c

clean:
	$(RM) *.o .nfs* core *.trace

clobber:
	$(RM) $(TARGETS)
