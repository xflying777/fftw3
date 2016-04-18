#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

void initial(fftw_complex *x, int N);
void print_complex_vector(fftw_complex *x, int N);

int main()
{
	printf("\n");

	int N = 8;
	fftw_complex *x;

	x = (fftw_complex *) malloc(N*sizeof(fftw_complex));
	printf(" fftw_complex mlloc success. \n");

	initial(x, N);
	print_complex_vector(x, N);

	return 0;
}

void initial(fftw_complex *x, int N)
{
	int i;

	for (int i=0; i<N; i++)
	{
		x[i][0] = 1.0*i;
		x[i][1] = 0.0;
	}

	printf(" Initial success. \n");
}

void print_complex_vector(fftw_complex *x, int N)
{
	int i;

	for (i=0; i<N; i++)
	{
		printf(" [%d] = %f , %f \n", i, x[i][0], x[i][1]);
	}
	printf("\n");
}

