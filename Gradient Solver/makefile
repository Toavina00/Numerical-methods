CC = g++

IFLAG = -Iinclude

DOBJ = obj
DINC = include
DSRC = src

_INC = Solver.hpp
_OBJ = main.o $(patsubst %.hpp, %.o, $(_INC))

OBJ = $(patsubst %, $(DOBJ)/%, $(_OBJ))
INC = $(patsubst %, $(DINC)/%, $(_INC))

$(DOBJ)/%.o: $(DSRC)/%.cpp $(INC)
	$(CC) -c $< -o $@ $(IFLAG)

solver: $(OBJ)
	$(CC) $^ -o $@


.PHONY: clean

clean:
	rm $(DOBJ)/*.o