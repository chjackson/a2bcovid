using namespace std;
#include "basicmodel.h"
#include "distributions.h"

void PreCalculateLikelihoods (run_params p, vector<double>& LNPreCalc, vector<double>& OGPreCalcP) {
	MakeLogNormalPreCalc (p,LNPreCalc);
	MakeOffsetGammaPreCalcP (p,OGPreCalcP);
}

void MakeLogNormalPreCalc (run_params p, vector<double>& LNPreCalc) {
	for (int i=0;i<500;i++) {
		//cout << i << " ";
		double L=LogNormal(i,p.smu,p.ssigma);
		//cout << L << "\n";
		LNPreCalc.push_back(L);
	}
}

double LogNormal (double x, double mu, double sigma) {
	double y=0;
	if (x<0) {
		y=-1e10;
	} else {
		if (x>=0.5) {
			y = gsl_cdf_lognormal_P(x+0.5,mu,sigma) - gsl_cdf_lognormal_P(x-0.5,mu,sigma);
		} else {
			y = gsl_cdf_lognormal_P(x+0.5,mu,sigma);
		}
		y=log(y);
	}
	return y;
}

double LogNormalP (double x, const vector<double>& LNPreCalc) {
	double y=-1e10;
	if (x>=0) {
		y=LNPreCalc[x];
	}
	return y;
}

void MakeOffsetGammaPreCalcP (run_params p, vector<double>& OGPreCalcP) {
	for (int i=-26;i<30;i++) {
		double L=OffsetGammaCDFFlex(i,p.pa,p.pb,p.po);
		OGPreCalcP.push_back(L);
	}
}

double OffsetGammaCDFFlex (double x, double a, double b, double o) {
	double y=0;
	if (x+o+0.5<0) {
		y=-1e10;
	} else {
		if (x+o-0.5<0) {
			y = gsl_cdf_gamma_P(x+o+0.5,a,b);
		} else {
			y = gsl_cdf_gamma_P(x+o+0.5,a,b)-gsl_cdf_gamma_P(x+o-0.5,a,b);
		}
		y=log(y);
	}
	return y;
}

double OffsetGammaPreCalcP (int x, const vector<double>& OGPreCalcP) {
	int z=x+26;
	double y=-1e10;
	if (z>=0&&z<OGPreCalcP.size()) {
		y=OGPreCalcP[z];
	}
	return y;
}

double Poisson (int n, double lambda) {
	double y = gsl_ran_poisson_pdf(n,lambda);
	y=log(y);
	return y;
}

void PrecalculateSymptomDist (run_params p, vector<int>& sdist_interval, vector<double>& sdist_prob) {
	for (int i=-3;i<30;i++) {
		sdist_interval.push_back(i);
		double y=exp(OffsetGammaCDFFlex (i, p.ucta, p.uctb, p.ucto));
		sdist_prob.push_back(y);
	}
}
