.PHONY: all compile exec clean
all: compile exec

compile:
	g++ -fopenmp ${opts} main${N}.cpp
exec:
	./a.out >/dev/null
clean:
	@rm a.out
	
