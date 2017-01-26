CXX=g++
CXXFLAGS=-Wall -c -g
LFLAGS=-Wall -g
SOURCES=main.cpp parser.cpp fit.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=CG

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $< -o $@