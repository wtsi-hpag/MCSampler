#pragma once
#include <iomanip>
#include <vector>
namespace MCMC
{

	struct ParameterEstimate
	{
		double Fraction;
		double Lower;
		double Upper;
		double Median;

		 friend std::ostream& operator<< (std::ostream& stream, const ParameterEstimate& p)
		 {
			stream << p.Median << "(" << p.Fraction << " CI = [-" << p.Lower << ", +" << p.Upper << "] )";
			return stream;
		 }
	};
	
	struct Surface
	{
		int Nx;
		int Ny;
		std::vector<double> X;
		std::vector<double> Y;
		std::vector<std::vector<double>> Z;

		Surface(int nx, int ny) : Nx(nx), Ny(ny)
		{
			X.resize(Nx);
			Y.resize(Ny);
			Z = std::vector<std::vector<double>>(Ny,std::vector<double>(Nx,0.0));
		}
	};

	struct Histogram
	{
		int Dimension;
		double LowerBound;
		double UpperBound;
		std::vector<double> RawData;
		std::vector<double> Centres;
		std::vector<double> Frequency;
		Histogram(){};
		Histogram(int bins)
		{
			Centres.resize(bins);
			Frequency.resize(bins,0.0);
		}
	};

}