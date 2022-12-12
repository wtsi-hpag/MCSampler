#pragma once
#include <vector>
#include <iostream>
#include <random>
#include <thread>
#include "acor.h"
#include "PredictionBar.h"
namespace MCMC
{
	
	const double NEG_INF = -99*1e300; 

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
		std::vector<WalkerEnsemble> Past;
		WalkerEnsemble Current;
		WalkerEnsemble PreviousState;
		std::vector<double> Means;
		int ThinningRate = 10;
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
	};
	struct Histogram
	{
		std::vector<double> Centres;
		std::vector<double> Frequency;
		Histogram(int bins)
		{
			Centres.resize(bins);
			Frequency.resize(bins,0.0);
		}
	};
	class Sampler
	{
		private:
			const int WalkerCount;
			const int Dimensions;
			const int ThreadCount;
			std::vector<std::string> ParameterNames;
			Walkers WalkerSet;
			std::default_random_engine generator;
			std::uniform_real_distribution<double> uniform = std::uniform_real_distribution<double>(0,1);
			std::uniform_int_distribution<int> pairSelector;

			bool Verbose = true;

			template<typename Functor>
			void SamplingStep(int kStart, int kEnd, Functor & LogScore)
			{
				std::vector<double> proposal(Dimensions);
				for (int k = kStart; k < kEnd; ++k)
				{
					int pair = pairSelector(generator);
					if (pair >= k)
					{
						pair += 1;
					}
					double zNum = uniform(generator) * (MoveParameter - 1) + 1;
					double Z = (zNum * zNum)/MoveParameter; //fancy inverse transform sampling on the GW10 function
					//samples from the function g(Z=z) \propto 1/sqrt(z) for z between 1/a and a

					for (int j = 0; j < Dimensions; ++j)
					{
						double pairJ = WalkerSet.Previous(pair,j);
						proposal[j] = pairJ + Z * (WalkerSet.Previous(k,j) - pairJ);
					}
					double newScore = LogScore(proposal);
					double oldScore = WalkerSet.PreviousScore(k);
					double acceptanceScore = std::min(1.0,exp((newScore - oldScore)));
					double r = uniform(generator);
					if (r<= acceptanceScore)
					{
						WalkerSet.Current.Positions[k] = proposal;
						WalkerSet.Current.Scores[k] = newScore;
					}
					
				}
			}
			void Comment(std::string input)
			{
				if (Verbose)
				{
					std::cout << input << std::endl;
				}
			}
		public:
			double MoveParameter = 2;
			Sampler(int nWalkers,int dimensions, int nThreads) : WalkerCount(nWalkers), Dimensions(dimensions), ThreadCount(nThreads)
			{
				Comment("Initialising Multi-Core Sampler, using default parameter names)");
				ParameterNames.resize(dimensions);
				for (int i = 0; i < dimensions; ++i)
				{
					ParameterNames[i] = "Param " + std::to_string(i);
				}
			}
			Sampler(int nWalkers, int dimensions) : WalkerCount(nWalkers), Dimensions(dimensions), ThreadCount(1)
			{
				Comment("Initialising Single-Core Sampler, using default parameter names)");
				ParameterNames.resize(dimensions);
				for (int i = 0; i < dimensions; ++i)
				{
					ParameterNames[i] = "Param " + std::to_string(i);
				}
			}
			Sampler(int nWalkers,int dimensions,int nThreads, std::vector<std::string> names) : WalkerCount(nWalkers), Dimensions(dimensions), ThreadCount(nThreads)
			{
				Comment("Initialising Multi-Core Sampler");
				int origSize = names.size(); //protect against dodgy name arrays
				if (origSize != dimensions)
				{
					Comment("\tWARNING: Name vector was ill-sized and had to be buffered");
					names.resize(dimensions);
					for (int j = origSize; j < dimensions; ++j)
					{
						names[j] = "Param " + std::to_string(j);
					}
				}
				ParameterNames = names;
			}
			Sampler(int nWalkers,int dimensions, std::vector<std::string> names) : WalkerCount(nWalkers), Dimensions(dimensions), ThreadCount(1)
			{
				Comment("Initialising Single-Core Sampler");
				int origSize = names.size(); //protect against dodgy name arrays
				if (origSize != dimensions)
				{
					Comment("\tWARNING: Name vector was ill-sized and had to be buffered");
					names.resize(dimensions);
					for (int j = origSize; j < dimensions; ++j)
					{
						names[j] = "Param " + std::to_string(j);
					}
				}
				ParameterNames = names;
			}
			template<typename Functor>
			void Run(Functor & f, int nSamples,std::vector<double> initialGuess)
			{
				Run(f,nSamples,1,initialGuess);
			}


			template<typename Functor>
			void Run(Functor & f, int nSamples,int thinningRate,std::vector<double> initialGuess)
			{
				if (initialGuess.size() != Dimensions)
				{
					initialGuess.resize(Dimensions,0.0);
				}

				WalkerSet = Walkers(nSamples,thinningRate,WalkerCount,Dimensions, initialGuess,generator);
				Comment("Walker Ensemble initialised\nPopulating initial scores....");		
				
				for (int k = 0; k < WalkerCount; ++k)
				{
					double score = f(WalkerSet.Current.Positions[k]);
					WalkerSet.Current.Scores[k] = score;
				}
				WalkerSet.Update();
				Comment("Initial conditions prepped, beginning main loop");

				pairSelector = std::uniform_int_distribution<int>(0,WalkerCount-2);//-2 to allow for the offset of not choosing yourself!
				int kPerThread = WalkerCount/ThreadCount;
				
				int nonMainCount = ThreadCount - 1;
				std::vector<std::thread> threads(nonMainCount);
				
				PredictionBar pb(nSamples);
				for (int l = 0; l < nSamples; ++l)
				{
					for (int t = 0; t < nonMainCount; ++t)
					{
						threads[t] = std::thread(&Sampler::SamplingStep<Functor>,this,t*kPerThread,(t+1)*kPerThread,std::ref(f));
					}
					SamplingStep(nonMainCount*kPerThread,WalkerCount,f);

					for (int t= 0; t < nonMainCount; ++t)
					{
						threads[t].join();
					}
					WalkerSet.Update();
					pb.Update(l);
				}
				Comment("Main loop complete, exiting normally");




			}

			void Seed(int n)
			{
				generator = std::default_random_engine(n);
			}

			Histogram GenerateHistogram(int burnIn, int dim, int bins)
			{
				Histogram out(bins);
				long int count = 0;
				std::vector<double> Series((WalkerSet.CurrentIdx - burnIn)*WalkerCount);

				for (int j = burnIn; j < WalkerSet.CurrentIdx; ++j)
				{
					for (int w = 0; w < WalkerCount; ++w)
					{
						double v = WalkerSet.Past[j].Positions[w][dim];
						Series[count] = v;
						++count;
					}
				}
				
				std::sort(Series.begin(),Series.end());
				// return out;
				double thresh = 0;
				int N = Series.size();
				int lowerIdx = thresh*N;
				int upperIdx = N - 1-lowerIdx;

				

			
				
				bool tailDominated = true;
				double lowestVal = Series[lowerIdx];
				double largestVal = Series[upperIdx];
				double delta;
				for (int its = 0; its < 4; ++its)
				{
					
					delta = (largestVal - lowestVal)/bins;

					// Histogram out(bins);
					for (int b = 0; b < bins; ++b)
					{
						out.Centres[b] = (b+0.5)*delta + lowestVal;
						out.Frequency[b] = 0;
					}


					for (int i =0; i < Series.size(); ++i)
					{
						double v= Series[i];
						if (v >= lowestVal && v <= largestVal)
						{
							int bin = (v - lowestVal)/delta;
							bin = std::min(bins-1,std::max(0,bin)); //needed because of equality in if statement: can cause under or overflows
							++out.Frequency[bin];
						}
					}

					double maxCont = 0;
					int bMax = -1;
					for (int b = 0; b < bins; ++b)
					{
						if (out.Frequency[b] > maxCont)
						{
							maxCont = out.Frequency[b];
							bMax = b;
						}
					}
					
					double contrast = 100;
					//trace up
					bool leftFound = false;
					bool rightFound = false;
					for (int b = 0; b < bMax; ++b)
					{
						double v = out.Frequency[b];
						double thisContrast = maxCont/(v + 1e-100);
						if (v > 0 && maxCont/v <= contrast)
						{
							lowestVal = out.Centres[b] - delta;
							b = bMax;
							leftFound = true;
						}
					}
					for (int b = bins-1; b > bMax; --b)
					{
						double v = out.Frequency[b];
						double thisContrast = maxCont/(v + 1e-100);
						if (v > 0 && maxCont/v <= contrast)
						{
							largestVal = out.Centres[b] + delta;
							b = bMax;
							rightFound = true;
						}
					}
					if (!leftFound)
					{
						lowestVal = out.Centres[bMax] - 2*delta;
					}
					if (!rightFound)
					{
						largestVal = out.Centres[bMax] + 2*delta;
					}
				}

				double r = 0;
				double prev = 0;
				for (int b = 0; b < bins; ++b)
				{
					out.Frequency[b] /= (count*delta);
					
					// r += out.Frequency[b];
					// out.Frequency[b] = r;
				}
				return out;
			}

			double ComputeAutocorrelation(int nWalkers)
			{
				int Duration = WalkerSet.CurrentIdx;

				double meanT = 0;

	

				for (int i = 0; i < nWalkers; ++i)
				{
					// 
					int T = WalkerSet.CurrentIdx;
					// std::cout << "Computing autocorr for " << i << " has " << T << " elements " << std::endl;
					std::vector<double> Series(T);

					double M = 0;
					for (int t = 0; t < T; ++t)
					{
						Series[t] = WalkerSet.Past[t].Scores[i];
						// M += Series[t];
					}
					double * copy = &Series[0];
					double Sigma;
					double Mean;
					double tau;
					int failure = acor(&Mean,&Sigma,&tau,copy,T);
					if (failure == 1)
					{
						std::cout << "Failed to accurately estimate autocor for chain of length " <<  T  << " on walker " << i << std::endl;
						tau = T/2;
					}
					tau = 1.0 + WalkerSet.ThinningRate * (tau -1);
					meanT += tau;
				}
				// std::cout << "Mean autocorr is " << meanT/nWalkers << std::endl;
				return meanT/nWalkers;
			}
	};
}