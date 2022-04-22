# dmft_jamming
Code for gradient descent in DMFT equations

Code used in https://arxiv.org/abs/2201.01161

Simulate the GD equations (5-10) in the paper.


Two methods:

  1) iterative solution (jamming_it):
		- starts with an initial guess for the kernels in Eqs. (6) k, MC, MR (or chi);
		- solve correlation-response Eqs. (8) at fixed kernels;
		- runs NP stochastic paths in Eqs. (5,7) at fixed kernels k, MC, chi and using Eqs. (8), for a range of h0 and time length N * dt;
		- evaluate the new kernels with Eq. (6) using the statistics on the stochastic paths.

  2) step-by-step solution (jamming_sbs):
		-	starts with the initial value of the kernels, analytic from Eq. (6);
		-	samples one time step in Eqs. (5,7,8) for a range of h0 and NP paths;
		-	compute the kernels from Eq. (6) using the statistics from the last step;
		-	repeat adding another time step until stop.


Compilation:
$METHOD = 'it' or 'sbs'


	gcc -c -Wall random.c
	gcc jamming_$METHOD.c random.o -Wall -o JAMMING_$METHOD -lm -lgsl -lgslcblas


Input parameters: see function [init_par](https://github.com/handrope/dmft_jamming/blob/d49b633c4e1b1fdd4845c8d729ecb30f0b96204e/jamming_it.c#L84) in each jamming_$METHOD.c
