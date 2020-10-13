#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include "Rcpp.h"
void PrintSeqDistances (const vector< vector<int> >& seqdists);
void PrintSequenceDistances (const vector< vector<int> >& seqdists);
void PrintVariants (const vector<sparseseq>& variants, const vector<pat>& pdat);
void PrintVariantLoci (const vector<allele>& allvar);
void PrintPdat (const vector<pat>& pdat);
void LikelihoodOutput (run_params p, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans);
DataFrame LikelihoodOutputR(run_params p, const vector<int>& ordered, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans);
void FinalOutput (run_params p, const vector<string>& removed, const vector<string>& fixed);
