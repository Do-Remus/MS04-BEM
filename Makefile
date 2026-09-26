# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -O2 -g -Wall -Wextra -std=c++17 -fno-omit-frame-pointer

# Output directory
OUT_DIR = out

# Target executable
TARGET = $(OUT_DIR)/helmholtz

# Source directory
SRCS_DIR = src

# Source files
SRCS = $(wildcard $(SRCS_DIR)/*.cpp)

# Main
MAIN = main.cpp

# Header directory
HEAD_DIR = src/headers

# Config directory
CONFIG_DIR = src/config

# Header files
HEADERS = $(wildcard $(HEAD_DIR)/*.hpp) $(wildcard $(CONFIG_DIR)/*.hpp)

# Objects directory
OBJ_DIR = out/obj

# Object files
OBJS = $(MAIN:%.cpp=$(OBJ_DIR)/%.o) $(SRCS:$(SRCS_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# ============================================================
# Profiling tools
# ============================================================

PERF = /usr/lib/linux-tools/5.4.0-216-generic/perf
FLAMEGRAPH = $(HOME)/FlameGraph

# Profiling output directory
PROFILE_DIR = profiling

# Perf / FlameGraph
PERF_DATA = $(PROFILE_DIR)/perf.data
PERF_SCRIPT = $(PROFILE_DIR)/perf.txt
PERF_FOLDED = $(PROFILE_DIR)/perf.folded
FLAMEGRAPH_SVG = $(PROFILE_DIR)/flamegraph.svg

# Memory profiling
MEMORY_REPORT = $(PROFILE_DIR)/memory.txt
CACHE_REPORT = $(PROFILE_DIR)/cache.txt
MASSIF_DATA = $(PROFILE_DIR)/massif.out
MASSIF_REPORT = $(PROFILE_DIR)/massif.txt
LEAK_REPORT = $(PROFILE_DIR)/leaks.txt


# ============================================================
# Default rule
# ============================================================

all: $(TARGET)


# ============================================================
# Build
# ============================================================

# Linking the object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compiling main.cpp
$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compiling source files
$(OBJ_DIR)/%.o: $(SRCS_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ============================================================
# Flame graph
# ============================================================

.PHONY: flamegraph

flamegraph: $(TARGET)
	@mkdir -p $(PROFILE_DIR)

	@echo "==> Recording CPU profile..."
	$(PERF) record -F 99 -g -o $(PERF_DATA) -- $(TARGET)

	@echo "==> Converting profile to text..."
	$(PERF) script -i $(PERF_DATA) > $(PERF_SCRIPT)

	@echo "==> Folding stack traces..."
	$(FLAMEGRAPH)/stackcollapse-perf.pl $(PERF_SCRIPT) > $(PERF_FOLDED)

	@echo "==> Generating flame graph..."
	$(FLAMEGRAPH)/flamegraph.pl $(PERF_FOLDED) > $(FLAMEGRAPH_SVG)

	@echo "==> Flame graph generated:"
	@echo "    $$(pwd)/$(FLAMEGRAPH_SVG)"

	@echo "==> Opening flame graph..."
	explorer.exe "$$(wslpath -w "$$(pwd)/$(FLAMEGRAPH_SVG)")" || true

# ============================================================
# Memory usage
# ============================================================

.PHONY: memory

memory: $(TARGET)
	@mkdir -p $(PROFILE_DIR)

	@echo "==> Measuring memory usage..."
	/usr/bin/time -v -o $(MEMORY_REPORT) $(TARGET)

	@echo "==> Memory report:"
	@cat $(MEMORY_REPORT)

# ============================================================
# Cache profiling
# ============================================================

.PHONY: cache

cache: $(TARGET)
	@mkdir -p $(PROFILE_DIR)

	@echo "==> Measuring cache behaviour..."

	$(PERF) stat \
		-o $(CACHE_REPORT) \
		-e cycles,instructions,\
cache-references,cache-misses,\
L1-dcache-loads,L1-dcache-load-misses,\
LLC-loads,LLC-load-misses,\
page-faults,minor-faults,major-faults \
		-- $(TARGET)

	@echo "==> Cache report:"
	@cat $(CACHE_REPORT)

# ============================================================
# Massif
# ============================================================

.PHONY: massif

massif: $(TARGET)
	@mkdir -p $(PROFILE_DIR)

	@echo "==> Running Massif..."

	valgrind \
		--tool=massif \
		--massif-out-file=$(MASSIF_DATA) \
		$(TARGET)

	@echo "==> Generating Massif report..."

	ms_print $(MASSIF_DATA) > $(MASSIF_REPORT)

	@echo "==> Massif report:"
	@cat $(MASSIF_REPORT)

# ============================================================
# Memory leaks
# ============================================================

.PHONY: leaks

leaks: $(TARGET)
	@mkdir -p $(PROFILE_DIR)

	@echo "==> Running Valgrind Memcheck..."

	valgrind \
		--leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		--log-file=$(LEAK_REPORT) \
		$(TARGET)

	@echo "==> Leak report:"
	@cat $(LEAK_REPORT)

# ============================================================
# Complete profiling
# ============================================================

.PHONY: profile

profile: flamegraph memory cache massif leaks
	@echo
	@echo "=========================================="
	@echo " Profiling terminé"
	@echo "=========================================="
	@echo "Résultats disponibles dans:"
	@echo "    $$(pwd)/$(PROFILE_DIR)/"

# ============================================================
# Clean
# ============================================================

clean:
	rm -f $(OBJS) $(TARGET)
	rm -rf $(PROFILE_DIR)

clean_output:
	rm -f ./outputs/*


# ============================================================
# Information
# ============================================================

show-headers:
	@echo "Fichiers d'en-tête :"
	@echo $(HEADERS)

show-objects:
	@echo "Fichiers objets :"
	@echo $(OBJS)

show-sources:
	@echo "Fichiers sources :"
	@echo $(SRCS)