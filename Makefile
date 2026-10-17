# Compiler
PLATFORM=$(shell uname -o)

ifeq ($(PLATFORM),Darwin)
#CXX = /opt/homebrew/Cellar/gcc/14.2.0/bin/g++-14 
CXX = /opt/homebrew/Cellar/gcc/16.2.0/bin/g++-16
CXXFLAGS = -Wc++11-extensions -std=c++14 -O3 -I /opt/homebrew/include/
else
CXX = g++
CXXFLAGS = -O6 -std=c++17 -I /home/aht/.cmdstan/cmdstan-2.35.0/stan/lib/stan_math/lib/boost_1.84.0 -fopenmp  
endif

# Source files
SRCS = test_arrestSeries.cpp poplr_arrest.cpp poplr.cpp read_utils.cpp 

# Object files
OBJS = $(SRCS:.cpp=.o)

all: poplr_arrest test_arrestSeries

poplr_arrest: poplr_arrest.o poplr.o read_utils.o arrest.o
	$(CXX) $(CXXFLAGS) -o poplr_arrest poplr_arrest.o poplr.o read_utils.o arrest.o

test_arrestSeries: test_arrestSeries.o poplr.o read_utils.o arrest.o
	$(CXX) $(CXXFLAGS) -o test_arrestSeries test_arrestSeries.o poplr.o read_utils.o arrest.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

poplr.o: poplr.hpp Makefile

read_utils.o: read_utils.hpp Makefile

poplr_arrest.o: arrestSeries.hpp read_utils.hpp Makefile read_utils.hpp arrest.hpp

test_arrestSeries.o: arrestSeries.hpp read_utils.hpp Makefile read_utils.hpp

clean:
	rm -f $(OBJS) poplr_arrest test_arrestSeries a.out
