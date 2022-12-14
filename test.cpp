#include "MCSampler.h"
#include "JSL.h"
double logLikelihood(const std::vector<double> & parameters)
{
	int dim = parameters.size();
	double r = 0;
	double sigma = 0.0001;
	for (int i = 0; i < dim; ++i)
	{
		double x = parameters[i];
		// if (x < -20)
		// {
		// 	return MCMC::NEG_INF;
		// }
		double d = (x - i*i)/(sigma);
		r +=   0.5*d*d;
	}
	return -r - 0.5 * dim * log(2*M_PI) - dim * log(sigma);
}

int main(int argc, char ** argv)
{
	int nWalkers = 70;//the number of members of the ensemble - the higher the number the more parallelisation helps you
	int dimensions = 20;//the number of parameters required for logLikelihodd
	int thin = 10;
	int threads = 1;
	if (argc > 1)
	{
		threads = std::stoi(argv[1]);
	}
	std::cout << "Running threads = " << threads << std::endl;
	MCMC::Sampler sampler(nWalkers, dimensions,threads);
	sampler.Seed(time(NULL));
	// sampler.BurnInFactor = 100;
	std::vector<double> init(dimensions,0.4);

	int nSamples = 1000000;
	sampler.BurnInFactor = 20;
	int worked = sampler.Run(logLikelihood, nSamples, thin,init);
	// if (worked == 1)
	// {
	// 	exit(2);
	// }
	
	double q = sampler.FunctionIntegrator();
	JSL::gnuplot gp;
	int plotDim = std::min(4,dimensions);
	gp.SetMultiplot(plotDim,1);
	gp.WindowSize(800,900);
	for (int d = 0; d < plotDim; ++d)
	{
		gp.SetAxis(d,0);
		int bins = 100;
		auto H = sampler.GenerateHistogram(d,bins);
		for (int b = 0; b < bins; ++b)
		{
			H.Frequency[b] *= pow(q,1.0/dimensions);
		}
		gp.Plot(H.Centres,H.Frequency);
		// gp.SetYLog(true);
		gp.SetGrid(true);
		
		std::cout << sampler.Estimate(d,0.69) << std::endl;
	}
	gp.Show();

	std::cout << pow(q,1.0/dimensions) << std::endl;
	// std::cout << "Function integration is: " << sampler.MeanFunction() << std::endl;

	return 0;
}