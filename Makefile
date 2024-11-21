CC	:= g++
INCDIRS	:=
LIBDIRS	:=
LIBS	:= -lcsdr++ -lfftw3f
CFLAGS	:= -O1 $(INCDIRS)

all: csdr-cwskimmer

csdr-cwskimmer: skimmer.o
	$(CC) $(CFLAGS) -o $@ $< $(LIBDIRS) $(LIBS)

clean:
	rm -f skimmer.o csdr-cwskimmer
