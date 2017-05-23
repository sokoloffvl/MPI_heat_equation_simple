// ConsoleApplication2.cpp:
//
//mpicc -o final final.c -lm
//mpirun -np 2 final
#include <stdio.h> 
#include <stdint.h>
#include <stdlib.h>
#include <math.h> 
#include <mpi.h>

int max_range;
int main(int argc, char* argv[])
{
	max_range = 100;
	   
	float M_2[max_range][max_range][max_range];
	float M_1[max_range][max_range][max_range];

	MPI_Init(&argc, &argv);
	MPI_Finalize();
	return 0;
}

