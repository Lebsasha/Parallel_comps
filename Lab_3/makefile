.PHONY: all exec compile hello

ifndef n
n=3
endif
ifndef opts
opts= -g -p -Wall
endif
all: compile exec
exec:
	mpirun -n ${n} bin/Debug/Temp
compile:
	mpic++ ${opts} main.cpp -o bin/Debug/Temp
hello:
	@echo ${CC}
