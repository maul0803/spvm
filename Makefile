# Compiler to use
CC = gcc

# Compiler flags
CFLAGS = -c -Wall

# Linker flags
LDFLAGS = -lopenblas

# Source files
SOURCES = timer.c spmv.c my_dense.c my_sparse.c

# Object files (change .c to .o)
OBJECTS = $(SOURCES:.c=.o)

# Executable name
EXECUTABLE = spmv

# Default target to build the executable
all: $(EXECUTABLE)

# Rule to link object files to create the final executable
$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS)

# Rule to compile each source file into an object file
%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

# Clean up generated files
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

# Rule to run the executable after building
run: $(EXECUTABLE)
	./$(EXECUTABLE)

.PHONY: all clean run

#make
#make run
#make clean
