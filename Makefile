
CC := cc

SRC := $(wildcard src/*.c)
OBJ := $(patsubst src/%.c,build/%.o,$(SRC))

CCFLAGS := -std=c11 -O3 -Wall -Wextra -Wpedantic -fopenmp
LDFLAGS := -lm -lSDL2 -fopenmp

TARGET := main

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $^ $(LDFLAGS) -o $@

build/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $< $(CCFLAGS) -c -o $@

clean:
	rm -rf $(TARGET) build

.PHONY: all clean

