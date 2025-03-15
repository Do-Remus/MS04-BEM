# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -std=c++11

# Output directory
OUT_DIR = out

# Target executable
TARGET = $(OUT_DIR)/sim202_program

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

# Objets directory
OBJ_DIR = out/obj

# Object files
OBJS = $(MAIN:%.cpp=$(OBJ_DIR)/%.o) $(SRCS:$(SRCS_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Default rule
all: $(TARGET)

# Linking the object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compiling source files into object files
$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -std=gnu++11
$(OBJ_DIR)/%.o: $(SRCS_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -std=gnu++11

# Clean rule to remove compiled files
clean:
	rm -f $(OBJS) $(TARGET)

show-headers:
	@echo "Fichiers d'en-tête :"
	@echo $(HEADERS)

show-objects:
	@echo "Fichiers objets :"
	@echo $(OBJS)

show-sources:
	@echo "Fichiers sources :"
	@echo $(SRCS)