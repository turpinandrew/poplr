# Compiler
CXX = /opt/homebrew/Cellar/gcc/14.2.0/bin/g++-14 

# Compiler flags
CXXFLAGS = -Wc++11-extensions -std=c++14 -O3 -I /opt/homebrew/include/

# Source files
SRCS = test_arrestSeries.cpp poplr_arrest.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

all: poplr test_arrestSeries

poplr: poplr_arrest.o
	$(CXX) $(CXXFLAGS) -o poplr poplr_arrest.o

test_arrestSeries: test_arrestSeries.o
	$(CXX) $(CXXFLAGS) -o test_arrestSeries test_arrestSeries.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

poplr_arrest.o: poplr_arrest.hpp arrestSeries.hpp

test_arrestSeries.o: poplr_arrest.hpp arrestSeries.hpp


clean:
	rm -f $(OBJS) poplr test_arrestSeries a.out