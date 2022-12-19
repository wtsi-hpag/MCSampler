#pragma once
#include <vector>
#include <random>
#include "acor.h"
#include "PredictionBar.h"
namespace MCMC
{
	double ale(double x, double y)
	{
		return std::max(x,y) + log(1.0 + exp(-abs(x - y)));
	}


	struct WalkerEnsemble
	{
		WalkerEnsemble(){}
		WalkerEnsemble(int walkers, int dimensions)
		{
			Positions = std::vector<std::vector<double>>(walkers, std::vector<double>(dimensions,0.0));
			Scores = std::vector<double>(walkers,0);
		}
		std::vector<std::vector<double>> Positions;
		std::vector<double> Scores;
	};


	struct Walkers
	{
		Walkers(){}
		Walkers(int duration, int thinningRate, int walkers, int dimension,std::vector<double> initialGuess,std::default_random_engine & generator)
		{
			WalkerCount = walkers;
			ThinningRate = thinningRate;
			Dimension = dimension;
			Past = std::vector<WalkerEnsemble>(duration/ThinningRate+1,WalkerEnsemble(walkers,dimension));
			CurrentTime = -1;
			Current = WalkerEnsemble(walkers,dimension);
			std::normal_distribution<double> updator(0,0.4);
			Means = std::vector<double>(walkers,-1);
		
			//generate initial state as a small gaussian ball around the single initial guess
			for (int i =0 ; i < dimension; ++i)
			{
				double a = initialGuess[i];
				if (a == 0)
				{
					a = 0.01;
				}
				for (int w = 0; w < walkers; ++w)
				{
					Current.Positions[w][i] = a * (1.0 + updator(generator));
				}
			}
			//don't update yet as scores not filled in yet!
		}

		void Expand(int newSamples)
		{
			
			Past.resize(Past.size() + newSamples/ThinningRate+1, WalkerEnsemble(WalkerCount,Dimension));
		}
		int CurrentTime;
		int CurrentIdx = 0;
		int WalkerCount;
		int Dimension;
		std::vector<WalkerEnsemble> Past;
		WalkerEnsemble Current;
		WalkerEnsemble PreviousState;
		std::vector<double> Means;
		int ThinningRate = 10;
		double ThinnedAutocorrelationTime;
		const WalkerEnsemble & Previous()
		{
			return PreviousState;
		}
		const double & Previous(int walker, int element)
		{
			return PreviousState.Positions[walker][element];
		}
		const double & PreviousScore(int walker)
		{
			return PreviousState.Scores[walker];
		}
		void Update()
		{
			++CurrentTime;
			// Past[CurrentTime] = Current;
			PreviousState = Current;
			if (CurrentTime % ThinningRate == 0)
			{
				Past[CurrentIdx] = Current;
				++CurrentIdx;
			}
		}
	

		double ComputeAutocorrelation()
		{
			int Duration = CurrentIdx;

			double meanT = 0;

			int failedAcors = 0;
			PredictionBar pb(WalkerCount);
			pb.SetName("\tProgress:");
			for (int i = 0; i < WalkerCount; ++i)
			{

				// 
				int T = CurrentIdx;
				// std::cout << "Computing autocorr for " << i << " has " << T << " elements " << std::endl;
				std::vector<double> Series(T);

				double M = 0;
				for (int t = 0; t < T; ++t)
				{
					Series[t] = Past[t].Scores[i];
				}
				double * copy = &Series[0];
				double Sigma;
				double Mean;
				double tau;
				int failure = acor(&Mean,&Sigma,&tau,copy,T);
				if (failure == 1 || std::isnan(T))
				{
					// std::cout << "Failed to accurately estimate autocor for chain of length " <<  T  << " on walker " << i << std::endl;
					tau = T;
					failedAcors++;
				}
				else
				{
					meanT += tau;
				}
				// std::cout << "\tWalker " << i << " reports " << tau << std::endl;
				// pb.Update(i);
			}
			
			double failureRate = (double)failedAcors/WalkerCount;
			std::cout << "Failure rate = " << failureRate << std::endl;
			double adjustedCorrelation;
			if (failureRate < 0.8)
			{
				ThinnedAutocorrelationTime = meanT/(WalkerCount - failedAcors);
				if (ThinnedAutocorrelationTime < 1.6)
				{
					ThinnedAutocorrelationTime = 1;
				}
				 adjustedCorrelation= std::max(ThinningRate * ThinnedAutocorrelationTime,1.0*ThinningRate);
			}
			else
			{
				ThinnedAutocorrelationTime = Duration;
				adjustedCorrelation = ThinningRate * Duration;
			}
			return adjustedCorrelation;
		}

		// double FunctionMean(double burnInFactor, int thinning)
		// {
		// 	double sum = 0;
		// 	int start = ThinnedAutocorrelationTime * burnInFactor;
		// 	long int c =0;
		// 	for (int i = start; i < CurrentIdx; i += thinning)
		// 	{
		// 		for (int w = 0; w < WalkerCount; ++w)
		// 		{
		// 			// sum += exp(Past[i].Scores[w]);
		// 			if (c == 0)
		// 			{
		// 				sum = Past[i].Scores[w];
		// 			}
		// 			else
		// 			{
		// 				sum = ale(sum,Past[i].Scores[w]);
		// 			}
		// 			++c;

		// 			double r = (double)rand()/RAND_MAX;

		// 		}
		// 	}
		// 	return exp(sum - log(c));
		// 	// return sum/c;
		// }
	};
}