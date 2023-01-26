CC = g++
CPPFLAGS=-std=c++14 -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor \
	  -Wold-style-cast -Wcast-align -Wuseless-cast \
	  -Wdouble-promotion -Wnull-dereference -Wmisleading-indentation \
	  -Wduplicated-cond -Wformat=2
DBGFLAGS = -g
OPTFLAG = -O3

BIN = simulation

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^

$(BIN): *.o
	$(CC) $(CPPFLAGS) -o $@ $^

debug: CPPFLAGS += $(DBGFLAGS)
debug : $(BIN)

release: CPPFLAGS += $(OPTFLAG)
release: $(BIN)

.PHONY: clean
clean:
	rm *.o
	rm $(BIN)
	rm *.dat
	rm *.gif

tidy:
	rm *.o

run:
	./$(BIN)

plot:
	./$(BIN)
	./*.plt
