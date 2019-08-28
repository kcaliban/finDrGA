# Specify compiler commands
CXX = g++ -O3 -Wall
# Specify name of file containing main function, without ext.
BINARY = PrisonersDilemma
# Get all library files, sorting removes duplicates
LIBS = $(sort $(basename $(wildcard lib*))) 
# Objects except binary and libraries
OBJ = $(basename $(filter-out $(addsuffix .cpp, $(LIBS)), $(wildcard *.cpp)))

all: compile

# Link everything together 
compile: libs objs
	$(CXX) $(addsuffix .o, $(OBJ)) -o $(BINARY) -L. $(addsuffix .so, $(LIBS))

# Generate shared libraries
libs: $(addsuffix .so, $(LIBS)) 

# Generate all object files except libraries
objs: $(addsuffix .o, $(OBJ))

# Generate shared library
%.so: %.cpp
	$(CXX) -fpic -shared -o libAtest.so $<

# Generate object file
%.o: %.cpp %.h
	$(CXX) -c $<

