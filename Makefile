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

PERF_DATA = perf.data
PERF_SCRIPT = perf.txt
PERF_FOLDED = perf.folded
FLAMEGRAPH_SVG = flamegraph.svg


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

flamegraph: $(TARGET)
	@echo "==> Recording CPU profile..."
	$(PERF) record -F 99 -g -- $(TARGET)

	@echo "==> Converting profile to text..."
	$(PERF) script > $(PERF_SCRIPT)

	@echo "==> Folding stack traces..."
	$(FLAMEGRAPH)/stackcollapse-perf.pl $(PERF_SCRIPT) > $(PERF_FOLDED)

	@echo "==> Generating flame graph..."
	$(FLAMEGRAPH)/flamegraph.pl $(PERF_FOLDED) > $(FLAMEGRAPH_SVG)

	@echo "==> Flame graph generated:"
	@echo "    $$(pwd)/$(FLAMEGRAPH_SVG)"

	@echo "==> Opening flame graph..."
	explorer.exe "$$(wslpath -w "$$(pwd)/flamegraph.svg")" || true


# ============================================================
# Clean
# ============================================================

clean:
	rm -f $(OBJS) $(TARGET)
	rm -f $(PERF_DATA) $(PERF_SCRIPT) $(PERF_FOLDED) $(FLAMEGRAPH_SVG)


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