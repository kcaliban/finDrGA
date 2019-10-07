# Specify compiler commands
CXX = mpic++ -O3 -Wall -fopenmp -std=c++11
# Specify name of file containing main function, without ext.
BINARY = Dvelopr
BINARY2 = PoolWorker
# Source files; ** wildcard does not work on my Make so just one level depth
SRCFILES = $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
# Filter out Test files
SRCFILES := $(filter-out $(wildcard src/*/*Test.cpp), $(filter-out $(wildcard src/*Test.cpp),$(SRCFILES)))
# Object files
# Change ext. and path
OBJFILES = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCFILES)))

# Testing
TEST = DveloprTest
SRCTEST = $(wildcard src/*Test.cpp) $(wildcard src/*/*Test.cpp)
OBJTEST = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCTEST)))

$(info    Source files: $(SRCFILES))
$(info    Object files: $(OBJFILES))
all: folders compile

# Create required folders
folders:
	mkdir -p obj
	mkdir -p obj/VinaInstance
	mkdir -p obj/PoolManager
	mkdir -p obj/GMXInstance
	mkdir -p obj/Serialization
# Link everything together 
compile: objs
	$(CXX) $(filter-out $(wildcard obj/PoolManager/PoolWorker.o),$(OBJFILES)) -o $(BINARY) -lm
	$(CXX) $(filter-out $(wildcard obj/Dvelopr.o),$(OBJFILES)) -o $(BINARY2) -lm

# Generate all object files
objs: $(OBJFILES)

# Test
test: testobjs
	$(CXX) $(OBJTEST) $(filter-out $(wildcard obj/PoolManager/PoolWorker.o), $(filter-out $(wildcard obj/Dvelopr.o),$(OBJFILES))) -o $(TEST) -lgtest -lgtest_main -lpthread
	./$(TEST) # --gtest_filter=Serialization.Pairs

testobjs: $(OBJTEST)

# Generate object file
obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $<

clean:
	rm -rf $(BINARY)
	rm -rf $(BINARY2)
	rm -rf obj
