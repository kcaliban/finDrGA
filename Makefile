# Specify compiler commands
CXX = g++ -O3 -Wall -fopenmp
# Specify name of file containing main function, without ext.
BINARY = PepGA
# Source files; ** wildcard does not work on my Make so just one level depth
SRCFILES = $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
# Filter out Test files
SRCFILES := $(filter-out $(wildcard src/*/*Test.cpp), $(filter-out $(wildcard src/*Test.cpp),$(SRCFILES)))
# Object files
# Change ext. and path
OBJFILES = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCFILES)))

# Testing
TEST = PepGATest
SRCTEST = $(wildcard src/*Test.cpp) $(wildcard src/*/*Test.cpp)
OBJTEST = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCTEST)))

$(info    Source files: $(SRCFILES))
$(info    Object files: $(OBJFILES))
all: folders compile test

# Create required folders
folders:
	mkdir -p obj
	mkdir -p obj/VinaInstance
	mkdir -p obj/PoolManager
	mkdir -p obj/GMXInstance
# Link everything together 
compile: objs
	$(CXX) $(OBJFILES) -o $(BINARY)

# Generate all object files
objs: $(OBJFILES)

# Test
test: testobjs
	$(CXX) $(OBJTEST) -o $(TEST) -lgtest -lgtest_main -lpthread
	./$(TEST)

testobjs: $(OBJTEST)

# Generate object file
obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $<

clean:
	rm -rf $(BINARY)
	rm -rf obj
