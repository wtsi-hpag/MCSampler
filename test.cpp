// #define GNUPLOT_NO_TIDY
#include "MCSampler.h"
#include "JSL.h"
double logLikelihood(const std::vector<double> & parameters)
{
	int dim = parameters.size();
	double r = 0;
	double r2 = 0;
	double sigma = 0.7;
	for (int i = 0; i < dim-1; ++i)
	{
		double x = parameters[i] + (i+1)*parameters[i+1]/300;
		// if (x < -20)
		// {
		// 	return MCMC::NEG_INF;
		// }
		double d = (x - i-2)/sigma;
		double d2 = (x + i)/sigma;
		r +=   0.5*d*d;
		r2 += 0.5 * d2*d2;
	}
	double x0 = (parameters[dim-1]);
	
	double a = -r;
	double b = log(1e-1) - r2;

	double q = std::max(a,b) + log(1.0 + exp(-abs(a-b)));
	return q - 0.5*x0*x0;
}

int main(int argc, char ** argv)
{
	int nWalkers = 40;//the number of members of the ensemble - the higher the number the more parallelisation helps you
	int dimensions = 3;//the number of parameters required for logLikelihodd
	int thin = 3;
	int threads = 1;
	if (argc > 1)
	{
		threads = std::stoi(argv[1]);
	}
	std::cout << "Running threads = " << threads << std::endl;
	MCMC::Sampler sampler(nWalkers, dimensions,threads);
	sampler.Seed(time(NULL));
	sampler.MoveParameter = 15;
	// sampler.BurnInFactor = 100;
	std::vector<double> init(dimensions,2.5);

	int nSamples = 1000000;
	// sampler.BurnInFactor = 0.3;
	int worked = sampler.Run(logLikelihood, nSamples, thin,init);
	// if (worked == 1)
	// {
	// 	exit(2);
	// }
	
	// auto S = sampler.GenerateSurface(0,1,100);
	JSL::gnuplot gp;
	int plotDim = std::min(4,dimensions);
	gp.SetMultiplot(plotDim,plotDim);
	gp.WindowSize(800,800);

	std::vector<MCMC::Histogram> Hs(plotDim);
	int bins = 100;
	for (int d = 0; d < plotDim; ++d)
	{
		
		
		auto H = sampler.GenerateHistogram(d,bins);
		std::cout << sampler.Estimate(d,0.69) << std::endl;
		Hs[d] = H;
	}
	gp.SetFontSize(JSL::Fonts::Global,7);
	// gp.SetFontSize(JSL::Fonts::Label,9);
	for (int d = 0; d < plotDim; ++d)
	{
		for (int j = 0; j < d; ++j)
		{
			gp.SetAxis(d,j);
			auto corr = sampler.GenerateCorrelationSurface(Hs[j],Hs[d],bins);
			gp.Map(corr.X,corr.Y,corr.Z);
			gp.SetXRange(Hs[j].LowerBound,Hs[j].UpperBound);
			gp.SetYRange(Hs[d].LowerBound,Hs[d].UpperBound);
			// gp.SetXLabel("")
		}
		gp.SetAxis(d,d);
		gp.Plot(Hs[d].Centres,Hs[d].Frequency);
		gp.SetXRange(Hs[d].LowerBound,Hs[d].UpperBound);
		// gp.SetYLog(true);
		gp.SetGrid(true);
		
	}
	gp.Show();

	// std::cout << "Function integration is: " << sampler.MeanFunction() << std::endl;

	return 0;
}