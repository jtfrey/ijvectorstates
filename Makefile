-include Makefile.inc

#

TARGET		= ftest

OBJECTS		= ijvectorstates.o ftest.o

MODULES		= ijvectorstates.mod

#

$(TARGET): $(OBJECTS)
	$(FC) -o $@ $+ $(LDFLAGS) $(LIBS)

clean:
	$(RM) $(MODULES) $(OBJECTS) $(TARGET)

ijvectorstates.o ijvectorstates.mod: ijvectorstates.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c ijvectorstates.F90

ftest.o: ijvectorstates.o ijvectorstates.mod ftest.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c ftest.F90

