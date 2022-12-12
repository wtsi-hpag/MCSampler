#pragma once
#include <vector>
#include <random>
#include "acor.h"
namespace MCMC
{
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
			Past = std::vector<WalkerEnsemble>(duration/ThinningRate+1,WalkerEnsemble(walkers,dimension));
			CurrentTime = -1;
			Current = WalkerEnsemble(walkers,dimension);
			std::normal_distribution<double> updator(0,0.01);
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
		int CurrentTime;
		int CurrentIdx = 0;
		int WalkerCount;
		
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
		double AutoCovariance(int walker, int delay)
		{
			
			if (Means[walker] == -1)
			{
				double sum3 = 0;
				for (int m = 0; m < CurrentIdx; ++m)
				{
					sum3 += Past[m].Scores[walker];
				}
				Means[walker] = sum3/CurrentIdx;
				std::cout << walker << " has mean f " << sum3/CurrentIdx << std::endl;
			}
			double mean = Means[walker];
			double sum = 0;
			
			for (int m = 0; m < CurrentIdx-delay; ++m)
			{
				double fLag = Past[m+delay].Scores[walker];
				double f = Past[m].Scores[walker];
				sum += (fLag - mean) * (f - mean);
			}
			return sum/(CurrentIdx - delay);
		}

		double ComputeAutocorrelation()
		{
			int Duration = CurrentIdx;

			double meanT = 0;



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
					// M += Series[t];
				}
				double * copy = &Series[0];
				double Sigma;
				double Mean;
				double tau;
				int failure = acor(&Mean,&Sigma,&tau,copy,T);
				if (failure == 1)
				{
					// std::cout << "Failed to accurately estimate autocor for chain of length " <<  T  << " on walker " << i << std::endl;
					tau = T;
				}
				tau = tau;
				meanT += tau;
			}
			ThinnedAutocorrelationTime = meanT/WalkerCount;
			double adjustedCorrelation = std::max(ThinningRate * ThinnedAutocorrelationTime,1.0*ThinningRate);
			return adjustedCorrelation;
		}
	};
}