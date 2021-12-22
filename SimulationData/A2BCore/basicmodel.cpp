#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>


#include "basicmodel.h"
#include "utilities.h"
#include "diagnostics.h"
#include "distributions.h"
#include "likelihoods.h"
#include "io.h"

int main(int argc, const char **argv) {

	int seed=(int) time(NULL);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);
	
	run_params p;
	GetOptions(p,argc,argv);
	
	vector<int> sdist_interval;
	vector<double> sdist_prob;
	PrecalculateSymptomDist(p,sdist_interval,sdist_prob);
	
	//Pre-calculate likelihoods
	vector<double> LNPreCalc; //LogNormal.  Assumes p.smu, p.ssigma parameters
	vector<double> OGPreCalcP; //Offset gamma
	PreCalculateLikelihoods(p,LNPreCalc,OGPreCalcP);

	int error=0;
	GetThresholds (p,p.threshold95,p.threshold99,p.t95NS,p.t99NS,error);
	if (error==1) {
		return 0;
	}

	//This code is set up just to run the key part of the likelihood code
	//Therefore does not read in CSV files etc; we know all of the parameters from a simulation of transmission
	//N.B. The code makes no use of location data
	
	//Input file: Use pat_file
	
	ifstream sim_file;
	sim_file.open(p.pat_file);
	int sa=0;
	int sb=0;
	int da=0;
	int db=0;
	int eva=0;
	int evb=0;
	for (int i=0;i<1000000;i++) {
		if (!(sim_file >> sa)) break;
		if (!(sim_file >> sb)) break;
		if (!(sim_file >> da)) break;
		if (!(sim_file >> db)) break;
		if (!(sim_file >> eva)) break;
		if (!(sim_file >> evb)) break;
		//Do calculation based on these values.
		
		//1.  Get locations for the pair.  Symptom onset dates are sa and sb
		//Range from min(sa,sb)-15 to max(sa,sb)+25
		vector<int> possible_tab;
		GetPossibleTimes (sa,sb,possible_tab);
		double lL=CalcLikelihood (p,sa,sb,da,db,eva,evb,possible_tab);
		//Here assume indiviudals in contact for half of the time.
		lL=log(exp(lL)/4.);
		cout << "Log L " << lL << " ";
		
		//Compare likelihood to threshold given D values
		ThresholdComparison (p,sa,sb,da,db,lL);
		
	}
	
	return 0;
}
	
