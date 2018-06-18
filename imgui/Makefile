CXX := g++
CC := gcc
LIB_OUT := imgui.a
CFLAGS := -Wall
#CXXFLAGS := -Wall -fPIC
CXXFLAGS := -Wall -O3

#C_OBJS := GL/gl3w.o
#    imgui_impl_opengl3.o
C_OBJS :=

CPP_OBJS := \
    imgui.o                 \
    imgui_demo.o            \
    imgui_draw.o            \
    imgui_impl_opengl2.o    \
    imgui_impl_sdl.o

INCLUDES :=  -I .

all: $(LIB_OUT) $(C_OBJS)

$(LIB_OUT): $(CPP_OBJS)
	$(AR) cr $(LIB_OUT) $(CPP_OBJS)

%.o:%.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) `sdl2-config --cflags` -c -o $@ $< `sdl2-config --libs` -lGL -ldl

GL/%.o:GL/%.c
	$(CC) $(INCLUDES) $(CFLAGS) `sdl2-config --cflags` -c -o $@ $< `sdl2-config --libs` -lGL -ldl

clean:
	@if [ -e imgui.o ]; then echo rm *.o; rm *.o; fi
	@if [ -e GL/gl3w.o ]; then echo rm GL/gl3w.o; rm GL/gl3w.o; fi
	@if [ -e imgui.a ]; then echo rm imgui.a; rm imgui.a; fi
	@if [ -e OUT ]; then echo rm OUT; rm OUT; fi
