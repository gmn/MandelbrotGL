
EXE := mandelbrot
CFLAGS_DEBUG := -g -Wall -std=c++11
CFLAGS_OPT := -Wall -std=c++11 -O3
INCLUDES := -Iimgui
LIBS :=
LINK := -lGL -lGLEW -lSDL2 -ldl
#LINK := -lGL -lSDL2 -ldl
IMGUI := imgui/imgui.a
#IMGUI := imgui/*.o
#GL3W := imgui/GL/gl3w.o
GL3W :=
CXX := g++
OBJS := mandelbrot.o

ifeq ($(MAKECMDGOALS),debug)
    override CFLAGS := $(CFLAGS_DEBUG)
else
    override CFLAGS := $(CFLAGS_OPT)
endif

.PHONY: debug

debug: $(EXE)

all: $(EXE)

mandelbrot: $(GL3W) $(OBJS) $(IMGUI)
	g++ $(CFLAGS) $(INCLUDES) $(LIBS) $^ $(LINK) -o $(EXE)

%.o:%.cpp
	$(CXX) $(INCLUDES) $(CFLAGS) `sdl2-config --cflags` -c -o $@ $< `sdl2-config --libs` -lGL -ldl

$(IMGUI):
	@cd imgui; make;

$(GL3W):
	@cd imgui; make;

clean:
	@if [ -e mandelbrot.o ]; then echo rm mandelbrot.o; rm mandelbrot.o; fi
	@if [ -e $(EXE) ]; then echo rm $(EXE); rm $(EXE); fi
	@if [ -e imgui.ini ]; then echo rm imgui.ini; rm imgui.ini; fi
	@if [ -e OUT ]; then echo rm OUT; rm OUT; fi

distclean: clean
	@cd imgui; make clean;

