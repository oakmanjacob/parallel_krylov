# Set basic parameters such as g++ compiler, 64 bit OS and compiler flags
BITS    ?= 64
CXX      = g++
LD       = g++
CXXFLAGS = -MMD -m$(BITS) -std=c++17 -Wall -Wextra -fPIC -I include -pg
LDFLAGS  = -m$(BITS) -lpthread -lspdlog -pg

# Set the out directory
ODIR      := ./build
OUTFOLDER := $(shell mkdir -p $(ODIR))

# All files with a main function
TARGETS = krylov

# All .cc files without a main function
CXXFILES = main gmres

# generate names of all .o and .exe files we're creaing
# name all .o files explicitly so can add to .PRECIOUS target
# .PRECIOUS prevents the .o files from being auto-removed
COMMONOFILES	= $(patsubst %, $(ODIR)/%.o, $(CXXFILES)) # linked to every executable
ALLOFILES 		= $(patsubst %, $(ODIR)/%.o, $(CXXFILES))
EXEFILES 		= $(patsubst %, $(ODIR)/%.exe, $(TARGETS))

# generate names of dependency files that g++ will generate, include later in makefile
DFILES = $(patsubst %.o, %.d, $(ALLOFILES))

.DEFAULT_GOAL = all
.PRECIOUS: $(ALLOFILES)
.PHONY: all debug clean

all: $(EXEFILES)
	@echo "Built for production"

# rules for building object files
$(ODIR)/%.o: src/%.cc
	@echo "[CXX] $< --> $@ (-O3)"
	@$(CXX) $< -o $@ -c $(CXXFLAGS) -O3

# rules for building executables
# assume an executable uses *all* of the common OFILES
$(ODIR)/%.exe: $(COMMONOFILES)
	@echo "[LD] $^ --> $@ (-O3)"
	@$(CXX) $^ -o $@ $(LDFLAGS) -O3 

# debug: $(EXEFILES)
# 	@echo "Built for debug"

# $(ODIR)/%.o: src/%.cc
# 	@echo "[CXX] $< --> $@ (-O0)"
# 	@$(CXX) $< -o $@ -c $(CXXFLAGS) -O0 -ggdb

# # rules for building executables
# # assume an executable uses *all* of the common OFILES
# $(ODIR)/%.exe: $(COMMONOFILES)
# 	@echo "[LD] $^ --> $@ (-O0)"
# 	@$(CXX) $^ -o $@ $(LDFLAGS) -O0

clean:
	@echo Cleaning up...
	@rm -rf $(ODIR)

# include any dependencies we generated previously
-include $(DFILES)