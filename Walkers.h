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
	struct WalkerScore
	{
		double Score=0;
		int ID;
	};

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
		Walkers(int duration, int thinningRate, int walkers, int dimension,std::vector<double> initialGuess,double confidence, std::default_random_engine & generator)
		{
			WalkerCount = walkers;
			ThinningRate = thinningRate;
			Dimension = dimension;
			Past = std::vector<WalkerEnsemble>(duration/ThinningRate+1,WalkerEnsemble(walkers,dimension));
			CurrentTime = -1;
			Current = WalkerEnsemble(walkers,dimension);
			double sigma = std::min(1000.0,std::max(1e-4,1.0/(confidence+1e-10)));
			std::normal_distribution<double> updator(0,sigma);

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
		int ViableCount;
		int Dimension;
		std::vector<WalkerEnsemble> Past;
		std::vector<bool> ChainViable;
		WalkerEnsemble Current;
		WalkerEnsemble PreviousState;
		// std::vector<double> Means;
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
	
		std::vector<WalkerScore> MeanScores()
		{
			std::vector<WalkerScore> Means(WalkerCount);
			int start = CurrentIdx/10;
			for (int w = 0; w < WalkerCount; ++w)
			{
				Means[w].ID = w;
				Means[w].Score = 0;
				for (int idx =start ; idx < CurrentIdx; ++idx)
				{
					Means[w].Score += Past[idx].Scores[w];
				}
			}
			for (int w = 0; w < WalkerCount; ++w)
			{
				Means[w].Score/=(CurrentIdx-start);
			}
			return Means;
		}

		double ComputeAutocorrelation()
		{
			int Duration = CurrentIdx;

			double meanT = 0;
			double meanAcor = 0;
			int failedAcors = 0;
			
			std::vector<double> taus;
			PredictionBar pb(WalkerCount);
			pb.SetName("\tProgress:");
			for (int i = 0; i < WalkerCount; ++i)
			{

				if (ChainViable[i])
				{
					int T = CurrentIdx;
					std::vector<double> Series(T);
					for (int t = 0; t < T; ++t)
					{
						Series[t] = Past[t].Scores[i];
						
					}
					
					
					double Sigma;
					double Mean;
					double tau;
					int failure = acor(&Mean,&Sigma,&tau,&Series[0],T);

					if (failure == 1 || std::isnan(tau))
					{
						if (!std::isnan(tau))
						{
							tau*=10;
						}
						else
						{
							tau = T;
						}
							// bool estimatorContinues = true;
							// double estimatorFactor = 7;
							// tau = 1.0;

							// double rho0 = 0;
							// for (int n = 0; n < T; ++n)
							// {
							// 	double x = Series[n] - Mean;
							// 	rho0 += x*x;
							// }
							// rho0/=T;
							// // std::cout << "For walker " << i << " rho = " << rho0 << std::endl;
							// double tauIdx = 1;
							// while (estimatorContinues)
							// {
							// 	double rho = 0;
							// 	for (int n = 0; n < T- tauIdx; ++n)
							// 	{
							// 		rho += (Series[n] - Mean) * (Series[n+tauIdx] - Mean);
							// 	}
							// 	rho/=(T - tauIdx);

							// 	tau += 2 * rho/rho0;
							
							// 	if (tauIdx >= estimatorFactor * tau || tauIdx >= T/2)
							// 	{
							// 		estimatorContinues = false;
							// 	}
							// 	++tauIdx;
							// }

					}


					taus.push_back(tau);
				}
				// std::cout << "\tWalker " << i << " reports " << tau << "  " << WalkerCount << "   \n" << std::endl;
				pb.Update(i+1);
			}

			if (taus.size() < 0.1 * WalkerCount)
			{
				ThinnedAutocorrelationTime = Duration;
				return Duration * ThinningRate;
			}
			std::sort(taus.begin(),taus.end());
			int middle = (taus.size() -1)/2;
			double medianTau = taus[middle];
			if (taus.size()%2 == 0)
			{
				medianTau = 0.5*medianTau + taus[middle+1];
			}
			

			double adjustedCorrelation;
			
			ThinnedAutocorrelationTime = medianTau;
			if (ThinnedAutocorrelationTime < 1.6)
			{
				ThinnedAutocorrelationTime = 1;
			}
			adjustedCorrelation= std::max(ThinningRate * ThinnedAutocorrelationTime,1.0*ThinningRate);
			
		
			return adjustedCorrelation;
		}

		void PruneChains()
		{
			std::vector<WalkerScore> means = MeanScores();
			std::sort(means.begin(),means.end(),[](const WalkerScore & lhs, const WalkerScore & rhs){return lhs.Score < rhs.Score;});
			double prev;
			double meanDiff = 0;
			for (int k = 1; k < means.size(); ++k)
			{
				meanDiff += means[k].Score - means[k-1].Score;
			}
			meanDiff /= (means.size()-1);


			std::vector<std::vector<int>> clusterMembers = {{means[0].ID}};
			std::vector<int> clusterBegin = {0};
			std::vector<int> clusterEnd;
			double fac = 10;
			int currentCluster = 0;
			for (int k = 1; k < means.size(); ++k)
			{
				double diff = (means[k].Score - means[k-1].Score);
				if (diff > fac * meanDiff)
				{
					clusterEnd.push_back(k);
					clusterBegin.push_back(k);
					clusterMembers.push_back({means[k].ID});
					++currentCluster;
				}
				else
				{
					clusterMembers[currentCluster].push_back(means[k].ID);
				}
			}
			if (clusterBegin.size()!=clusterEnd.size())
			{
				clusterEnd.push_back(means.size()-1);
			}
			// std::cout << " I have identified " << clusterBegin.size() << " clusters " << std::endl;
			// for (int j = 0; j < clusterMembers.size(); ++j)
			// {
			// 	std::cout << "\t" << clusterBegin[j] << " -> " << clusterEnd[j] << " members = {";
			// 	for (int q = 0; q < clusterMembers[j].size(); ++q)
			// 	{
			// 		clusterMembers
			// 	}
			// }
			ChainViable.resize(WalkerCount,true);
			ViableCount = WalkerCount;
			for (int j = 0; j < clusterMembers.size() - 1; ++j)
			{
				for (int q = 0; q < clusterMembers[j].size(); ++q)
				{
					int idx = clusterMembers[j][q];
					--ViableCount;
					ChainViable[idx] = false;
				}
			}
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