# Specify compiler commands
CXX = g++ -O3 -Wall -fopenmp
# Specify name of file containing main function, without ext.
BINARY = VinaGA
# Source files; ** wildcard does not work on my Make so just one level depth
SRCFILES = $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
# Object files
# Change ext. and path
OBJFILES = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCFILES)))

$(info    SRCFILES is $(SRCFILES))
$(info    OBJFILES is $(OBJFILES))
all: compile

# Link everything together 
compile: objs
	$(CXX) $(OBJFILES) -o $(BINARY)

# Generate all object files
objs: $(OBJFILES)

# Generate object file
obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $<

