
EXE := mandelbrot
CFLAGS := -g -Wall -std=c++11
CFLAGS_OPT := -Wall -std=c++11 -O3
INCLUDES := -ISDL2-2.0.8/include/
LIBS := -LSDL2-2.0.8/build/.libs
LINK := -lGL -lGLEW -lSDL2

all: $(EXE)

mandelbrot: mandelbrot.cpp
	g++ $(CFLAGS_OPT) $(INCLUDES) $(LIBS) mandelbrot.cpp $(LINK) -o $(EXE)

debug: mandelbrot.cpp
	g++ $(CFLAGS) $(INCLUDES) $(LIBS) mandelbrot.cpp $(LINK) -o $(EXE)

clean:
	@if [ -e $(EXE) ]; then echo rm $(EXE); rm $(EXE); fi

distclean: clean

