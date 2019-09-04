# Specify compiler commands
CXX = g++ -O3 -Wall -fopenmp
# Specify name of file containing main function, without ext.
BINARY = GMXTest
# Objects except binary and libraries
OBJ = $(basename $(wildcard *.cpp))

all: compile

# Link everything together 
compile: objs
	$(CXX) $(addsuffix .o, $(OBJ)) -o $(BINARY)

cleansub:
	rm work/*.top* work/*.edr* work/*.tpr* work/*.gro* work/*.itp* work/*.itp* work/*.log* work/*.trr* work/*.mdp* work/clean.pdb work/*.cpt*
# Generate all object files except libraries
objs: $(addsuffix .o, $(OBJ))

# Generate object file
%.o: %.cpp %.h
	$(CXX) -c $<

