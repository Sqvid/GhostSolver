CXX = g++
CXXFLAGS = -std=c++14 -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast \
	  -Wdouble-promotion -Wnull-dereference -Wmisleading-indentation \
	  -Wduplicated-cond -Wformat=2 -g -O0
LDFLAGS = -lm
OPTFLAG = -O2

BIN = main

OBJ = simulation.o eulerData.o slopeLimiter.o toroTests.o main.o

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^

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
