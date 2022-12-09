# MCSampler

MCSampler is an implementation of the 'Stretch-Move' algorithm for Markvo-Chain Monte Carlo sampling developed by Goodman & Weare (2010), written in a modern C++ template format for arbitrary functional analysis.

This work is heavily inspired by the Python implementation thereof, emcee, by Foreman-Mackey et. al (2013), and also makes use of the autocorrelation function Acor, also written by Goodman (https://www.math.nyu.edu/~goodman/software/acor/).

Parallelisation is supported through the ```<thread>``` library Unlike the Foreman-Mackey implementation, detailed balance is maintained on the entire ensemble by partioning the ensemble's history and its current position, allowing for arbitrary parellisation without impacting on the statistical validity of the ensemble.

## Usage

```c++
	MCMC::Sampler sampler(nWalkers, dimensions, nthreads)
	sampler.Run(targetFunction, nSamples, initialGuess)
```
<!-- 
## Compiling 

Currently only tested with C++17 specification on clang (Mac). Further testing needed. -->