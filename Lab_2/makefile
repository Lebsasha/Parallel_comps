.PHONY: all exec compile hello
all: compile exec
exec:
	mpirun -n ${n} bin/Debug/Temp
compile:
	mpic++ main.cpp -o bin/Debug/Temp
hello:
	@echo ${CC}
