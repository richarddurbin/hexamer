all: hexamer hextable

hexamer: hexamer.c readseq.c readseq.h
	cc -g -o hexamer hexamer.c readseq.c

hextable: hextable.c readseq.c readseq.h
	cc -g -o hextable hextable.c readseq.c -lm

clean:
	\rm  *.o hexamer hextable worm.hex *~
