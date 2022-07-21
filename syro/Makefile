CC=gcc
CFLAGS=-Wall
OBJS:=korg_syro_comp.o korg_syro_func.o korg_syro_volcasample.o korg_syro_volcasample_example.o
PROG:=korg_syro_volcasample_example

$(PROG): $(OBJS)
	$(CC) $^ -o $@ -lm

.PHONY: test clean

clean:
	-rm $(PROG) $(OBJS)

test: $(PROG)
	./$(PROG) out.wav s01:kick.wav
