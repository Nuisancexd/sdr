EXEC = sdr_app
CC = g++
CFLAGS = -Wall -O2 -std=c++17
LDFLAGS = -liio
SRC = main.cpp sdr.cpp

OBJ = $(SRC:.cpp=.o)
EXEC = laced

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)
