CC = gcc
CFLAGS = -Wall -O2 -std=c++17
LDFLAGS = -liio -lfftw3f -lm
SRC = main.cpp sdr.cpp FFT.cpp

OBJ = $(SRC:.cpp=.o)
EXEC = sdr

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)
