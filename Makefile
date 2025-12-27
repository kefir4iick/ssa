CC = gcc
CFLAGS = -Wall -Wextra -O2
TARGET = ring_test

all: $(TARGET)

$(TARGET): main.c ring.c
	$(CC) $(CFLAGS) -o $(TARGET) main.c ring.c

clean:
	rm -f $(TARGET) *.o

test: $(TARGET)
	./$(TARGET)
