#include "basicmodel.h"
#include <Rcpp.h>
#include <Rmath.h>

double OffsetGammaCDFFlex (double x, double a, double b, double o) {
	double y=0;
	if (x+o+0.5<0) {
		y=-1e10;
	} else {
		if (x+o-0.5<0) {
		  // y = gsl_cdf_gamma_P(x+o+0.5,a,b);
		  y = R::pgamma(x+o-0.5, a, b, true, false);
		} else {
//		  y = gsl_cdf_gamma_P(x+o+0.5,a,b)-gsl_cdf_gamma_P(x+o-0.5,a,b);
		  y = R::pgamma(x+o+0.5,a,b,true,false)  - R::pgamma(x+o-0.5,a,b,true,false);
		}
		y=log(y);
	}
	return y;
}

double LogNormal (double x, double mu, double sigma) {
	double y=0;
	if (x<0) {
		y=-1e10;
	} else {
		if (x>=0.5) {
//		  y = gsl_cdf_lognormal_P(x+0.5,mu,sigma) - gsl_cdf_lognormal_P(x-0.5,mu,sigma);
		  y = R::plnorm(x+0.5,mu,sigma, true, false) - R::plnorm(x-0.5,mu,sigma, true, false);
		} else {
//		  y = gsl_cdf_lognormal_P(x+0.5,mu,sigma);
		  y = R::plnorm(x+0.5,mu,sigma, true, false);
		}
		y=log(y);
	}
	return y;
}

double Poisson (int n, double lambda) {
//	double y = gsl_ran_poisson_pdf(n,lambda);
  double y = R::dpois(n,lambda,false);
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
