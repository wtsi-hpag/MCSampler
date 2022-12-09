# MCSampler

MCSampler is an implementation of the 'Stretch-Move' algorithm for Markvo-Chain Monte Carlo sampling developed by Goodman & Weare (2010), written in a modern C++ template format for arbitrary functional analysis.

This work is heavily inspired by the Python implementation thereof, emcee, by Foreman-Mackey et. al (2013), and also makes use of the autocorrelation function Acor, also written by Goodman (https://www.math.nyu.edu/~goodman/software/acor/).

Parallelisation is supported through the ```<thread>``` library Unlike the Foreman-Mackey implementation, detailed balance is maintained on the entire ensemble by partioning the ensemble's history and its current position, allowing for arbitrary parellisation without impacting on the statistical validity of the ensemble.

## Usage

The basic functionality of the code is as follows:
```c++
#include "MCSampler.h"
double logLikelihood(const std::vector<double> & parameters)
{
	return //some function of the input parameters
}

int main(int argc, char** argv)
{
	int nWalkers = //the number of members of the ensemble - the higher the number the more parallelisation helps you
	int dimensions = //the number of parameters required for logLikelihodd
	MCMC::Sampler sampler(nWalkers, dimensions;)
	sampler.Run(targetFunction, nSamples, initialGuess);
}
```

This works on any function which takes an ```std::vector<double>``` as input, and returns a double (including Lambdas). If the function requires more inputs -- such as a dataset which is being evaluated against -- then it is probably best to construct a functor which contains these parameters. If the functor possesses an overload for ```operator ()(std::vector<double>)```, then it works identically to the above:

```c++
#include "MCSampler.h"
class LikelihoodFunctor
{
	inputData MyData
	LiklihoodFunctor(inputData data)
	{
		MyData = data;
	}
	double operator ()(const std::vector<double> & parameters)
	{
		return //some function of the input parameters and MyData
	}
}
int main(int argc, char** argv)
{
	int nWalkers = //the number of members of the ensemble - the higher the number the more parallelisation helps you
	int dimensions = //the number of parameters required for logLikelihodd
	MCMC::Sampler sampler(nWalkers, dimensions;)

	inputData data = ///some data stored in some structure 

	LikelihoodFunctor targetFunction(data);
	sampler.Run(targetFunction, nSamples, initialGuess);
}

```
This produces an MCMC chain for each walker, computes the autocorrelation time, and then flattens the chains into a single, shuffled array of samples. Drawing samples from the distribution is then equivalent to taking subsequent numbers from this array.


### Things to do


- Estimate
- Plotter
- Draw samples
- Flattener
- Resume function?

<!-- 
## Compiling 

Currently only tested with C++17 specification on clang (Mac). Further testing needed. -->