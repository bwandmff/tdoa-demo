# Makefile for TDOA Algorithm Implementation
# High-Efficiency C Programming Standards

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -std=c99 -O2
LDFLAGS = -lm

# Directories
SRC_DIR = .
OBJ_DIR = obj
BIN_DIR = bin
TEST_DIR = .

# Source files
SOURCES = tdoa.c main.c
OBJECTS = $(SOURCES:%.c=$(OBJ_DIR)/%.o)

# Targets
TARGET = $(BIN_DIR)/tdoa_demo
TEST_TARGET = $(BIN_DIR)/tdoa_test

# Header dependencies
HEADERS = tdoa.h

# ============================================================================
# Build Rules
# ============================================================================

.PHONY: all clean test run benchmark help

all: $(TARGET)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR)/%.o: %.c $(HEADERS) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

# ============================================================================
# Testing
# ============================================================================

test: $(TARGET)
	@echo "Running TDOA demo..."
	./$(TARGET)

# Run with valgrind memory check
memtest: $(TARGET)
	@echo "Running memory check with valgrind..."
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes \
		--error-exitcode=1 ./$(TARGET)

# Run with clang static analysis
analyze:
	@echo "Running static analysis with clang-tidy..."
	clang-tidy -checks='*' tdoa.c main.c -- $(CFLAGS) -I.

# ============================================================================
# Development
# ============================================================================

run: test

# Debug build
debug: CFLAGS += -g -O0 -DDEBUG
debug: clean $(TARGET)

# Release build (optimized)
release: CFLAGS += -DNDEBUG -O3
release: clean all

# Benchmark
benchmark: $(TARGET)
	@echo "Running performance benchmark..."
	time ./$(TARGET)

# ============================================================================
# Utilities
# ============================================================================

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
	find . -name "*.o" -delete
	find . -name "*.d" -delete

# Count lines of code
loc:
	@echo "Lines of code:"
	wc -l tdoa.h tdoa.c main.c

# Format code
format:
	clang-format -i tdoa.c main.c

help:
	@echo "TDOA Algorithm Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build the demo (default)"
	@echo "  test      - Build and run demo"
	@echo "  memtest   - Run with valgrind memory check"
	@echo "  analyze   - Run clang-tidy static analysis"
	@echo "  debug     - Build with debug symbols"
	@echo "  release   - Build optimized release"
	@echo "  benchmark - Run performance benchmark"
	@echo "  clean     - Remove build artifacts"
	@echo "  loc       - Count lines of code"
	@echo "  format    - Format source code"
	@echo "  help      - Show this help message"
