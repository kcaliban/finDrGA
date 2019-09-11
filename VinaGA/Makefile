# Specify compiler commands
CXX = g++ -O3 -Wall -fopenmp -std=c++11
# Specify name of file containing main function, without ext.
BINARY = VinaGA
# Source files; ** wildcard does not work on my Make so just one level depth
SRCFILES = $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
SRCFILES := $(filter-out $(wildcard src/*/*Test.cpp), $(filter-out $(wildcard src/*Test.cpp),$(SRCFILES)))
# Object files
# Change ext. and path
OBJFILES = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCFILES)))

$(info    Source files: $(SRCFILES))
$(info    Object files: $(OBJFILES))
all: compile

# Link everything together 
compile: objs
	$(CXX) $(OBJFILES) -o $(BINARY)

# Generate all object files
objs: $(OBJFILES)

# Generate object file
obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $<

clean:
	rm $(BINARY)
	rm obj/*.o
	rm obj/*/*.o
