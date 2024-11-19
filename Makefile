CC	:= g++
INCDIRS	:=
LIBDIRS	:=
LIBS	:= -lcsdr++ -lfftw3f
CFLAGS	:= -O1 $(INCDIRS)

csdr-cwskimmer: skimmer.o
	$(CC) $(CFLAGS) -o $@ $< $(LIBDIRS) $(LIBS)

