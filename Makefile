CC = gcc
CFLAGS = -O0 -std=c++17
LDFLAGS = -liio -lfftw3f -lm -lstdc++ -lmatplot /home/clown/matplotplusplus/build/source/3rd_party/libnodesoup.a
SRC = main.cpp sdr.cpp FFT.cpp

OBJ = $(SRC:.cpp=.o)
EXEC = sdr

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)
