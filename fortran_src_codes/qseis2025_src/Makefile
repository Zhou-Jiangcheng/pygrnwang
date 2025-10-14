# Cross-platform Makefile for Windows and Linux
# Detect operating system
ifeq ($(OS),Windows_NT)
    DETECTED_OS := Windows
    EXE_EXT := .exe
    MKDIR := mkdir
    PATHSEP := \\
else
    DETECTED_OS := $(shell uname -s)
    EXE_EXT := .bin
    MKDIR := mkdir -p
    PATHSEP := /
endif

# Compiler settings
FC      := gfortran

# Base flags for both platforms
BASE_FFLAGS := -ffixed-line-length-none -std=legacy

# Platform-specific flags
ifeq ($(DETECTED_OS),Windows)
    FFLAGS := $(BASE_FFLAGS) -static -static-libgcc -static-libgfortran
else
    FFLAGS := $(BASE_FFLAGS) -fPIC -Wl,--no-relax
endif

# Source and object files
SRC     := $(wildcard *.f)
OBJ     := $(SRC:.f=.o)

# Target executables
TARGET_BIN := qseis2025$(EXE_EXT)
TARGET_DBG := qseis2025_dbg$(EXE_EXT)

.PHONY: all debug release clean info run

all: release

release: FFLAGS += -O3
release: $(TARGET_BIN)

debug: FFLAGS += -g -O0
debug: $(TARGET_DBG)

$(TARGET_BIN): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(TARGET_DBG): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

info:
	@echo "Detected OS: $(DETECTED_OS)"
	@echo "Compiler: $(FC)"
	@echo "Flags: $(FFLAGS)"
	@echo "Target: $(TARGET_BIN)"

run: $(TARGET_BIN)
ifeq ($(DETECTED_OS),Windows)
	$(TARGET_BIN)
else
	./$(TARGET_BIN)
endif

clean:
	rm -f *.o $(TARGET_BIN) $(TARGET_DBG)
