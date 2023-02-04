CXX = g++
CXXFLAGS = -std=c++14 -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast \
	  -Wdouble-promotion -Wnull-dereference -Wmisleading-indentation \
	  -Wduplicated-cond -Wformat=2
LDFLAGS = -lm
DBGFLAGS = -g
OPTFLAG = -O2

BIN = simulation

OBJ = fvmSimulation.o fvmData.o toroTests.o simulation.o

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^

debug: CXXFLAGS += $(DBGFLAGS)
debug : $(BIN)

release: CXXFLAGS += $(OPTFLAG)
release: $(BIN)

.PHONY: clean
clean:
	$(RM) *.o
	$(RM) $(BIN)
	$(RM) *.dat
	$(RM) *.gif

tidy:
	$(RM) *.o

run:
	./$(BIN)

plot:
	./*.plt
