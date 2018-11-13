PLAN_MODE=FFTW_PATIENT
CC=g++
CFLAGS=-DPLAN_MODE=$(PLAN_MODE) 
CFLAGS+=-c -Wall -std=c++11 -g -ggdb -fpermissive
LDFLAGS=-lfftw3
SOURCES= grueLattice.cpp \
						fft.cpp \
						circulant_ring.cpp \
						ksw_key.cpp \
						operations.cpp \
						rgsw.cpp \
						rlwe.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=grueLattice

all: $(SOURCES) $(EXECUTABLE)

check: 
	make check -C tests/

tests: 
	make -C tests/
	
stats: 
	make -C stats/

clean: 
	rm -f $(EXECUTABLE) $(OBJECTS)

clean_all: 
	rm -f $(EXECUTABLE) $(OBJECTS)
	make clean -C tests/
	make clean -C stats/
	rm -fr doc
	
doc: $(SOURCES)
	/usr/bin/doxygen Doxyfile
	make -C doc/latex

rebuild: clean;
	make -j5

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
