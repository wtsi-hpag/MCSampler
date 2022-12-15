#pragma once
#include <vector>
#include <iostream>
#include <random>
#include <thread>
#include <stdexcept>
#include "acor.h"
#include "PredictionBar.h"
#include "Walkers.h"

namespace MCMC
{
	
	const double NEG_INF = -99*1e300; 


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
			std::vector<bool> ThreadWaiting;
			std::vector<bool> OperationComplete;
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
					double acceptanceScore = std::min(1.0,pow(Z,Dimensions-1) * exp((newScore - oldScore)));
					double r = uniform(generator);
					if (r<= acceptanceScore)
					{
						WalkerSet.Current.Positions[k] = proposal;
						WalkerSet.Current.Scores[k] = newScore;
					}
					
				}
			}

			template<typename Functor>
			void ThreadWaiter(int threadID,int kStart, int kEnd, Functor & LogScore)
			{
				OperationComplete[threadID] = false;
				while (OperationComplete[threadID] == false)
				{
					if (ThreadWaiting[threadID])
					{
						SamplingStep(kStart,kEnd,LogScore);

						ThreadWaiting[threadID] = false;
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

			template<typename Functor>
			int MainSampleLoop(Functor & f, int nSamples)
			{
				Comment("\tBeginning main loop");
				pairSelector = std::uniform_int_distribution<int>(0,WalkerCount-2);//-2 to allow for the offset of not choosing yourself!
				int kPerThread = WalkerCount/ThreadCount;
				
				int nonMainCount = ThreadCount - 1;
				std::vector<std::thread> threads(nonMainCount);
				ThreadWaiting.resize(nonMainCount);
				OperationComplete.resize(nonMainCount);
				
				for (int i =0; i < nonMainCount; ++i)
				{
					threads[i] = std::thread(&Sampler::ThreadWaiter<Functor>,this,i,i*kPerThread,(i+1)*kPerThread,std::ref(f));
				}

				PredictionBar pb(nSamples);
				pb.SetName("\tProgress: ");
				for (int l = 0; l < nSamples; ++l)
				{
					for (int t = 0; t < nonMainCount; ++t)
					{
						//threads[t] = std::thread(&Sampler::SamplingStep<Functor>,this,t*kPerThread,(t+1)*kPerThread,std::ref(f));
						ThreadWaiting[t] = true;
					}
					SamplingStep(nonMainCount*kPerThread,WalkerCount,f);
					int t = 0;
					while (t < nonMainCount)
					{
						if (ThreadWaiting[t] == false)
						{
							++t;
						}
					}
					
					WalkerSet.Update();
					if (Verbose)
					{
						pb.Update(l);
					}
				}
				if (Verbose)
				{
					pb.Clear();
				}
				for (int t = 0; t < nonMainCount; ++t)
				{
					OperationComplete[t] = true;
					threads[t].join();
				}
				Comment("\tMain loop complete\n\tComputing Autocorrelation Time");

				double tau = WalkerSet.ComputeAutocorrelation();
				double burnIn = tau * BurnInFactor;

				if (burnIn < nSamples)
				{
					Comment("\tMean autocorrelation time found to be " + std::to_string(tau) + " with a pre-thinning rate of " +std::to_string(WalkerSet.ThinningRate));
					
					int newRate = ceil(  (WalkerSet.ThinnedAutocorrelationTime)-0.1);
					if (newRate > 1)
					{
						Comment("\t\tAdditional thinning by a factor of " + std::to_string(newRate) + " recommended");
						AdditionalThinningRate = newRate;
					}				
					
					Comment("\tSample density judged sufficient");
					return 0;
				}
				else
				{
					double rate = (double)tau/nSamples;
					Comment("\tMean autocorrelation time found to be " + std::to_string(tau));
					Comment("\n-----------------------------------------------------");
					Comment("\t\tWARNING!");
					Comment("\ttau = " + std::to_string(rate) + " x sampled positions.");
					Comment("\tSample not statistically valid"); //don't throw an error as can be resumed!
					Comment("-----------------------------------------------------");
					return 1;
				}
			}
		public:
			double MoveParameter = 2;
			double BurnInFactor = 5; //The number of autocorrelation times used as the burn in period. 
			int AdditionalThinningRate = 1;
			Sampler(int nWalkers,int dimensions, int nThreads) : WalkerCount(nWalkers), Dimensions(dimensions), ThreadCount(nThreads)
			{
				Comment("Initialising Multi-Core Sampler, using default parameter names");
				ParameterNames.resize(dimensions);
				for (int i = 0; i < dimensions; ++i)
				{
					ParameterNames[i] = "Param " + std::to_string(i);
				}
			}
			Sampler(int nWalkers, int dimensions) : WalkerCount(nWalkers), Dimensions(dimensions), ThreadCount(1)
			{
				Comment("Initialising Single-Core Sampler, using default parameter names");
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
			int Run(Functor & f, int nSamples,int thinningRate,std::vector<double> initialGuess)
			{
				Comment("\nNew MCMC Run Beginning");
				if (initialGuess.size() != Dimensions)
				{
					initialGuess.resize(Dimensions,0.0);
				}

				WalkerSet = Walkers(nSamples,thinningRate,WalkerCount,Dimensions, initialGuess,generator);
				Comment("\tWalker Ensemble initialised\n\tPopulating initial scores.");		
				
				for (int k = 0; k < WalkerCount; ++k)
				{
					double score = f(WalkerSet.Current.Positions[k]);
					WalkerSet.Current.Scores[k] = score;
				}
				WalkerSet.Update();
				

				
				return MainSampleLoop(f,nSamples);

			}

			template<typename Functor>
			int Resume(Functor & f, int nSamples)
			{
				Comment("Resuming previous operation with " + std::to_string(nSamples)+ " samples");

				WalkerSet.Expand(nSamples);
				Comment("\tWalker Set expanded");
				
				return MainSampleLoop(f,nSamples);

			}

			void Seed(int n)
			{
				generator = std::default_random_engine(n);
			}

			Histogram GenerateHistogram(int dim, int bins)
			{
				Histogram out(bins);
				long int count = 0;
				int thinBurn = WalkerSet.ThinnedAutocorrelationTime * BurnInFactor;
				std::vector<double> Series((WalkerSet.CurrentIdx - thinBurn)*WalkerCount/AdditionalThinningRate);

				// std::cout << "Computing hist for " << thinBurn << "  " << WalkerSet.CurrentIdx << std::endl;
				for (int j = thinBurn; j < WalkerSet.CurrentIdx; j+=AdditionalThinningRate)
				{
					for (int w = 0; w < WalkerCount; ++w)
					{
						double v = WalkerSet.Past[j].Positions[w][dim];
						Series[count] = v;
						++count;

					}
				}
				
				// std::sort(Series.begin(),Series.end());
				// return out;
				int N = Series.size();
				
				bool tailDominated = true;
				double lowestVal = *std::min_element(Series.begin(),Series.end());
				double largestVal = *std::max_element(Series.begin(),Series.end());
				double delta;
				for (int its = 0; its < 40; ++its)
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
					
					double contrast =100;
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
							// std::cout << "Left contrast = " << thisContrast << "  " << lowestVal << std::endl;
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



			std::vector<std::vector<double>> FlattenedChains(int thinningRate)
			{
				int burnIn = WalkerSet.ThinnedAutocorrelationTime * BurnInFactor;
				int size = (WalkerSet.CurrentIdx - burnIn)/thinningRate * WalkerCount;
				std::cout << "Flattened vector has size " << size << std::endl;
				std::vector<std::vector<double>> out(size,std::vector<double>(Dimensions));

				int c = 0;
				for (int i = burnIn; i < WalkerSet.CurrentIdx; i+=thinningRate)
				{
					for (int w = 0; w < WalkerCount; ++w)
					{
						out[c] = WalkerSet.Past[i].Positions[w];
						++c;
					}
				}

				return out;
			}
			std::vector<std::vector<double>> FlattenedChains()
			{
				return FlattenedChains(AdditionalThinningRate);
			}

			const std::vector<double> & DrawPosition()
			{
				int burnIn = WalkerSet.ThinnedAutocorrelationTime * BurnInFactor;
				int size = (WalkerSet.CurrentIdx - burnIn);

				int w = floor(uniform(generator) * WalkerCount);
				int t = burnIn + floor(uniform(generator) * size);

				return WalkerSet.Past[t].Positions[w];
			}
			std::vector<std::vector<double>> DrawPositions(int n)
			{
				int burnIn = WalkerSet.ThinnedAutocorrelationTime * BurnInFactor;
				int size = (WalkerSet.CurrentIdx - burnIn);
				if (n > size/AdditionalThinningRate)
				{
					throw std::invalid_argument("You requested more draws than there are independent samples in the chain, cannot comply with this request");
				}

				std::vector<std::vector<double>> out(n,std::vector<double>(Dimensions,0));

				for (int i =0; i < n; ++i)
				{
					int w = floor(uniform(generator) * WalkerCount);
					int t = burnIn + floor(uniform(generator) * size);
					out[i] = WalkerSet.Past[t].Positions[w];
				}
				return out;
			}

			ParameterEstimate Estimate(int dim, double fraction)
			{
				int thinBurn = WalkerSet.ThinnedAutocorrelationTime * BurnInFactor;
				std::vector<double> Series((WalkerSet.CurrentIdx - thinBurn)*WalkerCount/AdditionalThinningRate);
				int count = 0;
				// std::cout << "Computing hist for " << thinBurn << "  " << WalkerSet.CurrentIdx << std::endl;
				for (int j = thinBurn; j < WalkerSet.CurrentIdx; j+=AdditionalThinningRate)
				{
					for (int w = 0; w < WalkerCount; ++w)
					{
						double v = WalkerSet.Past[j].Positions[w][dim];
						Series[count] = v;
						++count;
					}
				}
				std::sort(Series.begin(),Series.end());
				int N = count;

				double upFrac = 0.5 + fraction/2;
				double midFrac = 0.5;
				double downFrac = 0.5 - fraction/2;
				std::vector<double> fracs = {downFrac, midFrac, upFrac};
				for (int i = 0; i < fracs.size(); ++i)
				{
					int low = floor(fracs[i] * N);
					double bruch = fracs[i] * N - low;
					double val = Series[low] + bruch * (Series[low+1] - Series[low]);
					fracs[i] = val;
				}
				ParameterEstimate p;
				p.Fraction = fraction;
				p.Median = fracs[1];
				p.Lower = p.Median - fracs[0];
				p.Upper = fracs[2] - p.Median;
				

				return p;
			}

			double FunctionIntegrator()
			{

				// return WalkerSet.FunctionMean(BurnInFactor,AdditionalThinningRate);
				std::vector<double> SemiMajor(Dimensions,0.0);
				std::vector<double> Centre(Dimensions,0.0);

				double nby2 = (double)Dimensions/2;
				double logA  = -nby2 * log(2*M_PI);
				for (int d = 0; d < Dimensions; ++d)
				{
					ParameterEstimate p = Estimate(d,0.5);
					Centre[d] = p.Median;
					SemiMajor[d] = std::min(abs(p.Lower),abs(p.Upper));
					std::cout << "\t" << d << "  " << Centre[d] << "  " << SemiMajor[d] << std::endl;
					logA -= log(SemiMajor[d]);
				}

				


				int thinBurn = WalkerSet.ThinnedAutocorrelationTime * BurnInFactor;
				double run = -9999999;
				int c = 0;
				for (int j = thinBurn; j < WalkerSet.CurrentIdx; j+= 1)
				{
					for (int w = 0; w < WalkerCount; ++w)
					{
						double r = 0;
						for (int d = 0; d < Dimensions; ++d)
						{
							double x = WalkerSet.Past[j].Positions[w][d];
							double q = (x - Centre[d])/SemiMajor[d];
							r -= 0.5 * q * q;
						}

					
						run = ale(run,logA + r- WalkerSet.Past[j].Scores[w]);
					
						++c;
					}
				}

				return  exp(log(c)-run);
			}

	};
}