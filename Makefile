# Compiler to use
CC = gcc

# Default build mode
BUILD = release

# Compiler flags for release mode
ifeq ($(BUILD),release)
    CFLAGS = -c -Wall -Wextra #-Ofast #Ofast for  speed (less accuracy)
else
    CFLAGS = -g -O0 -c -Wall -Wextra #g for Debug
endif

# Linker flags
LDFLAGS = -lopenblas

# Libraries
# LDLIBS = -lm -lgslcblas -lgsl
LDLIBS = -lm -lopenblas -lgsl

# Source files
SOURCES = timer.c spmv.c my_dense.c my_sparse.c

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Executable name
EXECUTABLE = spmv

# Default target to build the executable
all: $(EXECUTABLE)

# Rule to link object files to create the final executable
$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

# Implicit rule to compile each source file into an object file
# %.o: %.c
# 	$(CC) $(CFLAGS) -c $< -o $@

# Clean up generated files
clean:
	#rm -f $(OBJECTS) $(EXECUTABLE)
	$(RM) $(OBJECTS) $(EXECUTABLE)

# Rule to run the executable after building
run: $(EXECUTABLE)
	./$(EXECUTABLE)

.PHONY: all clean run # targets aren't files

# make
# make run
# make clean
# make BUILD=debug
# make BUILD=release