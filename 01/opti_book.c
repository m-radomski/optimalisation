#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <assert.h>

#define FALSE 0
#define TRUE 1

#include "real_problem.h"

uint32_t state = 0;

uint32_t xorshift32()
{
	uint32_t x = state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	state = x;
	return x;
}

double uniform_random(double a, double b) 
{
	uint32_t rand_int = xorshift32();
	rand_int = xorshift32();
	double r = (double)(rand_int) / (double)(0xffffffff);

	double arg_span = b - a;
    return a + (r * arg_span);
}

double target_func(double x)
{
	return -cos(0.1 * x) * exp(-(pow((0.1*x - 2*M_PI), 2))) + 0.002 * pow((0.1*x), 2);
}

double target_func2(double x)
{
	return fabs(sym_get_max_temp(x) - 50.0);
}

typedef struct
{
	double start_x0;
	double start_x1;

	double expansion_coeff;
	double a;
	double b;

	int success;
} ExpansionResult;

ExpansionResult expand_region(double b, double (*tfunc)(double), double expansion_coeff)
{
	int N_MAX = 10000;
	double x0 = 0.0;
	double x1 = uniform_random(x0, b);

	ExpansionResult result = { 0 };
	result.start_x0 = x0;
	result.start_x1 = x1; 
	result.expansion_coeff = expansion_coeff;
	result.success = TRUE;

	if(tfunc(x0) == tfunc(x1))
	{
		// printf("First out: |%05.05f|%05.05f|\n", x0, x1);
		result.a = x0;
		result.b = x1;
		return result;
	}

	if(tfunc(x1) > tfunc(x0))
	{
		x1 = -x1;
		if(tfunc(x1) >= tfunc(x0))
		{
			// printf("Second out: |%05.05f|%05.05f|\n", x0, x1);
			result.a = x0;
			result.b = x1;
			return result;
		}
	}

	int i = 0;
	double *xs = malloc(N_MAX * 2 * sizeof(double));
	xs[0] = x0;
	xs[1] = x1;
	do
	{
		if(i > N_MAX)
		{
			result.success = FALSE;
			return result;
		}

		i += 1;
		xs[i+1] = pow(expansion_coeff, i) * x1;

	} while(!(tfunc(xs[i]) <= tfunc(xs[i+1])));

	if(tfunc(xs[i-1]) < tfunc(xs[i+1])) {
		result.a = xs[i-1];
		result.b = xs[i+1];
		return result;
	} else {
		result.a = xs[i+1];
		result.b = xs[i-1];
		return result;
	}

	assert(FALSE);
}

double minimize_fibonacci(double a, double b, double accuracy, double (*tf)(double))
{
	double phi[32] = { 0 };
	int fibonacci_number_num = 32;

	double L = (b-a) / accuracy;
	int k = 0;
	
	for(int i = 0; i < fibonacci_number_num; i++)
	{
		if(i == 0) {
			phi[i] = 0.0;
		} else if(i == 1) {
			phi[i] = 1.0;
		} else {
			phi[i] = phi[i-1] + phi[i-2];
		}
	}

	for(int i = 0; i < fibonacci_number_num; i++) {
		if(phi[i] > L) { k = i; break; }
	}

	double *A = malloc(k * sizeof(double));
	double *B = malloc(k * sizeof(double));
	double *C = malloc(k * sizeof(double));
	double *D = malloc(k * sizeof(double));

	A[0] = a;
	B[0] = b;
	C[0] = B[0] - (phi[k-1]/phi[k])*(B[0] - A[0]);
	D[0] = A[0] + B[0] - C[0];

	for(int i = 0; i <= k - 4; i++)
	{
		if(tf(C[i]) < tf(D[i]))
		{
			// Constrain to the left
			A[i+1] = A[i];
			B[i+1] = D[i];
		}
		else
		{
			// Constrain to the right
			B[i+1] = B[i];
			A[i+1] = C[i];
		}

		C[i+1] = B[i+1] - (phi[k-i-2]/phi[k-i-1])*(B[i+1] - A[i+1]);
		D[i+1] = A[i+1] + B[i+1] - C[i+1];
	}

	double res = C[k - 3];

	free(A);
	free(B);
	free(C);
	free(D);

	return res; // x*
}

double square(double x)
{
	return x*x;
}

double minimize_lagrange(double a, double c, double b, double epsilon_accuracy, double gamma_accuracy, double (*tf)(double))
{
	int n = 100;

	double *A = malloc((n+5) * sizeof(double));
	double *B = malloc((n+5) * sizeof(double));
	double *C = malloc((n+5) * sizeof(double));
	double *D = malloc((n+5) * sizeof(double));

	A[0] = a;
	B[0] = b;
	C[0] = c;

	int i = 0;
	for(i = 0; i < n; i++)
	{
		double numerator =
			tf(A[i])*(square(C[i])-square(B[i])) +
			tf(C[i])*(square(B[i])-square(A[i])) +
			tf(B[i])*(square(A[i])-square(C[i]));
		double denominator =
			tf(A[i])*(C[i]-B[i]) +
			tf(C[i])*(B[i]-A[i]) +
			tf(B[i])*(A[i]-C[i]);

		D[i] = 0.5 * numerator / denominator;

		if(A[i] < D[i] < C[i])
		{
			if(tf(D[i]) < tf(C[i]))
			{
				A[i + 1] = A[i];
				C[i + 1] = D[i];
				B[i + 1] = C[i];
			}
			else
			{
				A[i + 1] = D[i];
				C[i + 1] = C[i];
				B[i + 1] = B[i];
			}
		}
		else
		{
			if(C[i] < D[i] < B[i])
			{
				if(tf(D[i]) < tf(C[i]))
				{
					A[i + 1] = C[i];
					C[i + 1] = D[i];
					B[i + 1] = B[i];
				}
				else
				{
					A[i + 1] = A[i];
					C[i + 1] = C[i];
					B[i + 1] = D[i];
				}
			}
			else
			{
				goto ERROR_RETURN_INFINITY;
			}
		}
	}

	if(i == n) {
		goto ERROR_RETURN_INFINITY;
	}

	double res = D[i];
	goto FREE_ALL;

ERROR_RETURN_INFINITY:
	res = INFINITY;

FREE_ALL:
	free(A);
	free(B);
	free(C);
	free(D);

	return res; // x*
}

int main(int argc, char **argv)
{
	/* 
	struct timeval t = { 0 };
	gettimeofday(&t, NULL);
	state = (t.tv_usec + t.tv_sec);

	int human_readable = argc > 1 && argv[1][0] == 'h';

	double a = -100.0;
	double b = 100.0;
	double expansion_coeff = 1.5;
	for(int i = 0; i < 10000; i++)
	{
		ExpansionResult r = expand_region(b, target_func, expansion_coeff);

		if(r.success){
			double fmin = minimize_fibonacci(r.a, r.b, 0.1, target_func);
			double lmin = minimize_lagrange(r.a, r.b, (r.b-r.a) * 0.5, 0.1, 0.1, target_func);

			char *fmt = human_readable ? 
				"|% 8.3lf|% 8.3lf|% 8.3lf|% 8.3lf|% 8.3lf|% 8.3lf|\n" :
				"%lf,%lf,%lf,%lf,%lf,%lf\n";
			printf(fmt, r.start_x0, r.start_x1, r.a, r.b, fmin, lmin);
		}
	}
	*/

	printf("%lf\n", sym_get_max_temp(24.127701));

	int human_readable = argc > 1 && argv[1][0] == 'h';
	double a = 1.0;
	double b = 100.0;
	double expansion_coeff = 1.5;
	for(int i = 0; i < 10; i++)
	{
		//ExpansionResult r = expand_region(b, target_func2, expansion_coeff);

		double fmin = minimize_fibonacci(1.0, 100.0, 0.001, target_func2);
		double lmin = minimize_lagrange(1.0, 100.0, (100.0-1.0) * 0.5, 0.1, 0.1, target_func2);

		char *fmt = human_readable ? 
			"|% 8.3lf|% 8.3lf|% 8.3lf|% 8.3lf|% 8.3lf|% 8.3lf|\n" :
			"%lf,%lf,%lf,%lf,%lf,%lf\n";
		printf(fmt, 1.0, 100.0, 1.0, 100.0, fmin, lmin);
	}
}
