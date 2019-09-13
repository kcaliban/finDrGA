# Specify compiler commands
CXX = g++ -O3 -Wall -fopenmp
# Specify name of file containing main function, without ext.
BINARY = PepGA
# Source files; ** wildcard does not work on my Make so just one level depth
SRCFILES = $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
SRCFILES := $(filter-out $(wildcard src/*/*Test.cpp), $(filter-out $(wildcard src/*Test.cpp),$(SRCFILES)))
# Object files
# Change ext. and path
OBJFILES = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCFILES)))

$(info    Source files: $(SRCFILES))
$(info    Object files: $(OBJFILES))
all: clean folders compile

# Create required folders
folders:
	mkdir obj
	mkdir obj/VinaInstance
	mkdir obj/PoolManager
	mkdir obj/GMXInstance
# Link everything together 
compile: objs
	$(CXX) $(OBJFILES) -o $(BINARY)

# Generate all object files
objs: $(OBJFILES)

# Generate object file
obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $<

clean:
	rm -rf $(BINARY)
	rm -rf obj
