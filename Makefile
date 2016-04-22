CXX = g++
CXXFLAGS = -O3 -Wall
CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

all: ngsParalog

.PHONY: clean

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp > $*.d

ngsParalog: $(OBJ)
	$(CXX) $(CXXFLAGS) -o ngsParalog *.o

clean:
	rm -f ngsParalog *.o *.d
