# Makefile for DGD-Cas9's connection_to_matrix binary (pipeline Step 6).
#
# NOTE: connection_to_matrix.cpp is NOT present in this branch's history --
# it needs to be supplied (from an earlier backup, or rewritten) before
# `make` will succeed. This Makefile just restores the documented build
# command from the README; it does not invent the missing C++ source.

CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17
TARGET   := connection_to_matrix
SRC      := connection_to_matrix.cpp

.PHONY: all clean check-source

all: check-source $(TARGET)

check-source:
	@if [ ! -f $(SRC) ]; then \
		echo "ERROR: $(SRC) not found."; \
		echo "This source file is missing from the repository history --"; \
		echo "see the README's 'Missing files' note. Supply it here before running 'make'."; \
		exit 1; \
	fi

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
