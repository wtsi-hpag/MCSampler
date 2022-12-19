// #define GNUPLOT_NO_TIDY
#include "JSL.h"
#include "MCSampler.h"

double logLikelihood(const std::vector<double> & parameters)
{
	int dim = parameters.size();
	double r = 0;
	double r2 = 0;
	double sigma = 3;
	for (int i = 0; i < dim-1; ++i)
	{
		double x = parameters[i] + parameters[i+1]*0.9;
		// if (x < -20)
		// {
		// 	return MCMC::NEG_INF;
		// }
		double d = (x - 6)/sigma;
		double d2 = (x + 3*i)/sigma;
		r +=   0.5*d*d;
		r2 += 0.5 * d2*d2;
	}
	double x0 = (parameters[dim-1])/sigma;
	
	double a = -r;
	double b = log(2e-1) - r2;

	double q = std::max(a,b) + log(1.0 + exp(-abs(a-b)));
	return q - 0.5*x0*x0;
}

int main(int argc, char ** argv)
{
	int nWalkers = 200;//the number of members of the ensemble - the higher the number the more parallelisation helps you
	int dimensions = 5;//the number of parameters required for logLikelihodd
	int thin = 20;
	int threads = 1;
	if (argc > 1)
	{
		threads = std::stoi(argv[1]);
	}
	std::cout << "Running threads = " << threads << std::endl;
	MCMC::Sampler sampler(nWalkers, dimensions,threads);
	sampler.Seed(time(NULL));
	sampler.MoveParameter = 400;
	// sampler.BurnInFactor = 100;
	std::vector<double> init(dimensions,-3);

	int nSamples = 1000000;
	sampler.StartingConfidence = 0.01;
	sampler.BurnInFactor = 50;
	int worked = sampler.Run(logLikelihood, nSamples, thin,init);
	// if (worked == 1)
	// {
	// 	exit(2);
	// }
	
	// auto S = sampler.GenerateSurface(0,1,100);
	

	
	int bins = 100;
	

	auto gp = sampler.CornerPlot(bins);
	gp.SetFontSize(JSL::Fonts::Global,7);
	// gp.SetFontSize(JSL::Fonts::Label,9);
	
	gp.Show();

	// std::cout << "Function integration is: " << sampler.MeanFunction() << std::endl;

	return 0;
}