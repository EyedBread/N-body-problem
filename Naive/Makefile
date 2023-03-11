all: Sequential Parallel

Sequential: SequentialN.cpp
	g++ SequentialN.cpp -o Sequential

Parallel: Npar.cpp
	g++ Npar.cpp -o Parallel -fopenmp

clean:
	rm -f Sequential Parallel