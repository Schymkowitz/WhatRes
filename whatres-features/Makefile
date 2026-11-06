# Makefile for whatres-features, the WhatRes structure-quality commands
# Switch Lab 2026
#
#   make            optimised binary ./whatres-features (default, native architecture)
#   make debug      binary with debug symbols ./whatres-features_debug
#   make clean      remove build artefacts
#   make install    copy ./whatres-features to PREFIX (default ~/bin)
#
# macOS architecture (optional): make ARCH=arm64 | ARCH=x86_64 | ARCH=universal
#
# Compiler: clang++ on macOS, g++ on Linux (any g++ with C++14 support).
# Every .cpp file in this directory is compiled.

TARGET    := whatres-features
BUILD_DIR := build

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    CXX := clang++
else
    CXX := g++
endif

ARCH ?=
ARCH_FLAGS :=
ifeq ($(UNAME), Darwin)
    ifeq ($(ARCH), universal)
        ARCH_FLAGS := -arch arm64 -arch x86_64
    else ifneq ($(strip $(ARCH)),)
        ARCH_FLAGS := -arch $(ARCH)
    endif
endif

SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

CXXFLAGS_COMMON  := -std=gnu++14 -I. -Wall -Wno-unused-parameter -fpermissive $(ARCH_FLAGS)
CXXFLAGS_RELEASE := $(CXXFLAGS_COMMON) -O3
CXXFLAGS_DEBUG   := $(CXXFLAGS_COMMON) -O0 -g

# Compiled atom-parameter library (see atomparams.h): lib/macos (arm64 + x86_64)
# or lib/linux-<arch>.
ifeq ($(UNAME), Darwin)
    AP_LIB := lib/macos/libatomparams.a
else
    AP_LIB := lib/linux-$(shell uname -m)/libatomparams.a
endif

# pdb_quality_check and structure_survey run on several threads
ifeq ($(UNAME), Linux)
    LDFLAGS := $(AP_LIB) -lpthread
else
    LDFLAGS := $(AP_LIB) $(ARCH_FLAGS)
endif

.PHONY: all
all: CXXFLAGS := $(CXXFLAGS_RELEASE)
all: $(TARGET)

.PHONY: debug
debug: CXXFLAGS := $(CXXFLAGS_DEBUG)
debug: $(TARGET)_debug

$(TARGET)_debug: $(OBJS) $(AP_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

$(TARGET): $(OBJS) $(AP_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)
	@echo "Built: $@  ($(UNAME) / $(CXX))"

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(AP_LIB):
	@echo "missing $(AP_LIB): no compiled parameter library for this platform"; exit 1

-include $(DEPS)

PREFIX ?= $(HOME)/bin

.PHONY: install
install: all
	install -d $(PREFIX)
	install -m 755 $(TARGET) $(PREFIX)/$(TARGET)
	@echo "Installed to $(PREFIX)/$(TARGET)"

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET) $(TARGET)_debug
