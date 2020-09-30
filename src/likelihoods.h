#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

#include "Rcpp.h"

struct placetime {
	int date;
	string ward;
};

bool comparetprob(tprob v1, tprob v2);

void CalculateTDLikelihoods (run_params p, const vector<pat>& pdat, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, vector< vector<ijlike> >& like_trans);
void RemoveDuplicateTimes (vector<tprob>& contact_times_probs);
void FillTimes (vector<tprob>& contact_times_probs);



double LikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, const vector<pat>& pdat);
double NoSeqLikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector<pat>& pdat);





void FindLikelihoodSubsetsIJ(run_params p, const vector< vector<ijlike> >& like_trans, vector< vector<int> >& subsets);


Rcpp::String ThresholdsR (run_params p, double L);
Rcpp::String ThresholdsNSR (run_params p, double L);

void Thresholds (run_params p, double L);
void ThresholdsNS (run_params p, double L);

void SetThreshold (run_params& p);

