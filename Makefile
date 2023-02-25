CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast \
	  -Wdouble-promotion -Wnull-dereference -Wmisleading-indentation \
	  -Wduplicated-cond -Wformat=2 -g -O0
LDFLAGS = -lm
OPTFLAG = -O3

BIN = main

OBJ = simulation.o eulerData.o slopeLimiter.o toro2DTests.o main.o

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^

release: CXXFLAGS += $(OPTFLAG)
release: $(BIN)

DOCDIRS = html/ latex/

.PHONY: clean
clean:
	$(RM) *.o
	$(RM) $(BIN)
	$(RM) *.dat
	$(RM) *.gif
	$(RM) -r $(DOCDIRS)

tidy:
	$(RM) *.o

rebuild: clean
rebuild: $(BIN)

run:
	./$(BIN)

plot:
	./$(BIN)
	./*.plt 2> /dev/null

docs:
	doxygen Doxyfile
