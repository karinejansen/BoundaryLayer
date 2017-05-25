/*Program to predict the boundary layer flow in a 2-d channel using an integral boundary layer method and conservations of mass (assuming a Pohlhausen profile)*/
#pragma warning (disable : 4996)
//#define _CRT_SECURE_NO_DEPRECATE
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <time.h>
//#include <errno.h>

#define M_PI 3.14159265358979323846264338327

int main()
{
	float rho, mu, U, h1, h2, L, eps, nu, *x, alfa, beta, gamma, *U0, *delta2, delta2old, *tau_w, *delta1, error, *lambda;
	int n, i, N;


	rho = 1.19;
	mu = 1.82e-5;  //dynamic viscosity
	U = 0.5;
	h1 = 1.2;
	h2 = 0.77;
	L = 4.0;
	eps = 10e-6;
	nu = mu / rho;  //kinematic viscosity

	/*
	rho = 1.161;
	nu = 15.89e-6;  //kinematic viscosity
	mu = nu * rho;  //dynamic viscosity
	U = 25;
	h1 = 5;
	h2 = 5;
	L = 0.5;
	eps = 10e-6;
	*/

	printf("input number of discretisations along x - n [-]\n");
	scanf_s("%d", &n);
	N = n + 1;

	//Open a pipe to send input to gnuplot
	FILE* gp = _popen("gnuplot -persistent", "w");
	// set title and axes titles
	fprintf(gp, "set title 'n=%d, h1=%.2f and h2=%.2f'\n", n, h1, h2);
	fprintf(gp, "set title font ', 20'\n");
	fprintf(gp, "set ylabel 'boundary layer thickness [m]'\n");
	fprintf(gp, "set xlabel 'location x [m]'\n");
	//Start a plot sequence with gnuplot '-' command
	fprintf(gp, "plot '-' with lines title 'delta1' \n");

	//start the clock
	clock_t start_t, end_t;
	double total_t;
	start_t = clock();
	printf("Starting of the program, start_t = %ld\n", start_t);

	printf("Input:\n rho = %10.4f [kg/m^3]\n mu = %10.8f [kg/(s m)]\n U = %10.4f [m/s]\n h1 = %10.4f [m]\n h2 = %10.4f [m]\n L = %10.4f [m]\n eps = %10.8f [-]\n n = %d [-]\n", rho, mu, U, h1, h2, L, eps, n);

	x = malloc(N*sizeof(float));
	U0 = malloc(N*sizeof(float));
	delta2 = malloc(N*sizeof(float));
	tau_w = malloc(N*sizeof(float));
	delta1 = malloc(N*sizeof(float));
	lambda = malloc(N*sizeof(float));

	x[0] = 0;			// intitial values
	U0[0] = U;
	delta2[0] = 0;
	lambda[0] = 0;
	//tau_w[0] = +INFINITY;
	delta1[0] = 0;
	alfa = (3. / 10.) - (lambda[0] / 120.);
	beta = (37. / 315.) - (lambda[0] / 945.) - (lambda[0] * lambda[0] / 9072.);
	gamma = 2. + (lambda[0] / 6.);
	tau_w[0] = beta*gamma*((mu*U0[0]) / delta2[0]);

	// iteration over every n
	for (i = 1; i < N; ++i)
	{
		x[i] = i * L / n;
		//x[i] = L / 2 * cos(i * M_PI / n + M_PI) + L / 2;
		delta2[i] = delta2[i - 1];

		do          // iteration for U0, lambda and delta2 at position x[i]
		{
			delta2old = delta2[i];
			U0[i] = (h1*U) / (h1 - (x[i] * (h1 - h2) / L) - ((alfa * delta2[i]) / beta));
			lambda[i] = ((delta2[i] * delta2[i])*(U0[i] - U0[i - 1])) / ((x[i] - x[i - 1])*nu*(beta*beta));
			alfa = (3. / 10.) - (lambda[i] / 120.);
			beta = (37. / 315.) - (lambda[i] / 945.) - (lambda[i] * lambda[i] / 9072.);
			gamma = 2. + (lambda[i] / 6.);
			delta2[i] = sqrt(((delta2[i - 1] * delta2[i - 1]) + 2. * nu*((beta * gamma) / U0[i])*(x[i] - x[i - 1])) / (1. + 2. * (2. + alfa / beta)*(U0[i] - U0[i - 1]) / U0[i]));
			error = (delta2[i] - delta2old) / delta2[i];

		} while (sqrt(error*error) > eps); //stop iteration when delta2old and delta2 new converged to convergence tolerance (epsilon)

		tau_w[i] = beta*gamma*((mu*U0[i]) / delta2[i]);			//calc tau_w and delta1 at x[i] with just found U0, lambda and delta2.
		delta1[i] = alfa*delta2[i] / beta;

	}
	//printf("xi\n");

	FILE *f;
	f = fopen("result.txt", "w");
	//errno_t fopen_s(FILE *f "result.txt", "w");
	//f = fopen_s(f, "result.txt", "w");

	for (i = 0; i < N; ++i)	// print all values and send each data point to GNUPLOT
	{
		//Send a new data point to GNUPLOT
		fprintf(gp, "%g %g \n", x[i], delta1[i]);
		// print in command window
		//printf("i = %2.2d xi = %2.4f U0(x) = %6.6e lambda(x) = %6.6e delta1(x) = %6.6e tau_w(x) = %6.6e \n", i, x[i], U0[i], lambda[i], delta1[i], tau_w[i]);
		// print to text file
		printf("%i\t%f\t%f\t%f\t%f\t%f\n", i, x[i], U0[i], lambda[i], delta1[i], tau_w[i]);
		//fprintf(f, "%2.2d %2.4f %6.6e %6.6e %6.6e %6.6e \n", i, x[i], U0[i], lambda[i], delta1[i], tau_w[i]);
	}
	for (i = 1; i < N; ++i)	// print all values and send each data point to GNUPLOT
	{
		fprintf(f, "%i\t%f\t%f\t%f\t%f\t%f\n", i, x[i], U0[i], lambda[i], delta1[i], tau_w[i]);
		//printf("%f\t%f\t%f\t%f\t%f\t%f\n", i, x[i], U0[i], lambda[i], delta1[i], tau_w[i]);
		//printf("%2.4f\n", x[i]);
	}
	fclose(f);



	/*
	printf("delta1(x)\n");
	for (i = 0; i < N; ++i)	// to copy to excell
	{
	printf("%6.6e\n", delta1[i]);
	} */

	//Send the e command indicating end of data point list
	fprintf(gp, "e");
	//Flush the pipe, i.e. send all the commands to gnuplot
	fflush(gp);

	free(x);
	free(U0);
	free(delta2);
	free(tau_w);
	free(delta1);
	//free(k);
	free(lambda);

	end_t = clock();
	printf("End of the program, end_t = %ld\n", end_t);
	total_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
	printf("Total time taken by CPU: %f\n", total_t);

	//Wait until user presses any key (second time, after entering n first) to terminate
	//getchar();
	//getchar();

	//testPrint1();

	return 0;
}

