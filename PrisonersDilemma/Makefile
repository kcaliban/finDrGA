# Specify compiler commands
CXX = g++ -O3 -Wall
# Specify name of file containing main function, without ext.
BINARY = PrisonersDilemma
# Objects except binary and libraries
OBJ = $(basename $(wildcard *.cpp))

all: compile

# Link everything together 
compile: objs
	$(CXX) $(addsuffix .o, $(OBJ)) -o $(BINARY)

# Generate all object files except libraries
objs: $(addsuffix .o, $(OBJ))

# Generate object file
%.o: %.cpp %.h
	$(CXX) -c $<

