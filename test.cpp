#include "MCSampler.h"
#include "JSL.h"
double logLikelihood(const std::vector<double> & parameters)
{
	int dim = parameters.size();
	double r = 0;
	for (int i = 0; i < dim; ++i)
	{
		double d = (parameters[i]-i);
		r +=   0.5*d*d;
	}
	return -r;
}

int main(int argc, char ** argv)
{
	int nWalkers = 700;//the number of members of the ensemble - the higher the number the more parallelisation helps you
	int dimensions = 6;//the number of parameters required for logLikelihodd
	int thin = 100;
	int threads = 3;
	if (argc > 1)
	{
		threads = std::stoi(argv[1]);
	}
	std::cout << "Running threads = " << threads << std::endl;
	MCMC::Sampler sampler(nWalkers, dimensions,threads);
	
	// sampler.BurnInFactor = 100;
	std::vector<double> init(dimensions,0);

	int nSamples = 10000;
	int worked = sampler.Run(logLikelihood, nSamples, thin,init);
	// if (worked == 1)
	// {
	// 	exit(2);
	// }
	
	// JSL::gnuplot gp;
	// gp.SetMultiplot(dimensions,1);
	// gp.WindowSize(800,900);
	// for (int d = 0; d < dimensions; ++d)
	// {
	// 	gp.SetAxis(d,0);
	// 	int bins = 100;
	// 	auto H = sampler.GenerateHistogram(d,bins);
	// 	gp.Plot(H.Centres,H.Frequency);
	// 	// gp.SetYLog(true);
	// 	gp.SetGrid(true);
		
	// }


	worked = sampler.Resume(logLikelihood,nSamples * 10);
	if (worked == 1)
	{
		exit(2);
	}
	
	// for (int d = 0; d < dimensions; ++d)
	// {
	// 	gp.SetAxis(d,0);
	// 	int bins = 100;
	// 	auto H = sampler.GenerateHistogram(d,bins);
	// 	gp.Plot(H.Centres,H.Frequency);
	// 	// gp.SetYLog(true);
	// 	gp.SetGrid(true);
		
	// }
	// gp.Show();
	return 0;
}