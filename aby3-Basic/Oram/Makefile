CC = g++
CFLAGS = -Wall -g -std=c++14
SRCDIR = ./src
BINDIR = ./bin
lIBFILE = ./oram.a

SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(SRCS))

all: clean lib

main: $(OBJS)  
	$(CC) $(CFLAGS) $(OBJS) -o $(BINDIR)/main

lib: main
	ar rvs $(lIBFILE) $(BINDIR)/*.o
  
$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BINDIR)/*.o $(BINDIR)/main ./*.a

run:
	./bin/main