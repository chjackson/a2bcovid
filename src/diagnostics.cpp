#include "basicmodel.h"
#include "likelihoods.h"

#include "Rcpp.h"

using namespace std;
using namespace Rcpp;

void PrintSeqDistances (const vector< vector<int> >& seqdists) {
	Rcout << "Sequence distances\n";
	for (int i=0;i<seqdists.size();i++) {
		for (int j=0;j<seqdists[i].size();j++) {
			Rcout << seqdists[i][j] << " ";
		}
		Rcout << "\n";
	}
}

void PrintSequenceDistances (const vector< vector<int> >& seqdists) {
	Rcout << "Sequence distance matrices\n";
	for (int i=0;i<seqdists.size();i++) {
		for (int j=0;j<seqdists[i].size();j++) {
			Rcout << seqdists[i][j] << " ";
		}
		Rcout << "\n";
	}
	Rcout << "\n";
}

void PrintVariants (const vector<sparseseq>& variants, const vector<pat>& pdat) {
	Rcout << "Variants\n";
	for (int i=0;i<variants.size();i++) {
		Rcout << i << " " << pdat[i].code << " " << pdat[i].time_s << " ";
		for (int j=0;j<variants[i].locus.size();j++) {
			Rcout << variants[i].locus[j] << variants[i].allele[j] << " ";
		}
		Rcout << "\n";
	}
}

void PrintVariantLoci (const vector<allele>& allvar) {
	Rcout << "Variant loci\n";
	for (int i=0;i<allvar.size();i++) {
		Rcout << allvar[i].loc << allvar[i].nuc << " ";
	}
	Rcout << "\n";
}

void PrintPdat (const vector<pat>& pdat) {
	for (int i=0;i<pdat.size();i++) {
		Rcout << pdat[i].code << " " << pdat[i].time_s << " ";
		for (int j=0;j<pdat[i].locat.size();j++) {
			Rcout << pdat[i].locat[j].ward << " " << pdat[i].locat[j].date << " " << pdat[i].locat[j].prob << " ";
		}
		Rcout << "\n";
	}
}

DataFrame LikelihoodOutputR(run_params p, const vector<int>& ordered, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans) {
    // count the number of entries in the output, so we can preallocate memory
    int k, n;
    k = 0;
    for (int i=0;i<like_trans.size();i++) {
    	for (int j=0;j<like_trans[i].size();j++) {
	      ++k;
	    }
    }
    n = k;
    CharacterVector from(n);
    CharacterVector to(n);
    NumericVector from_hcw(n);
    NumericVector to_hcw(n);
    NumericVector orderi(n);
    NumericVector orderj(n);
	NumericVector likelihood(n);
    CharacterVector consistency(n);
    LogicalVector under_threshold(n);
    k = 0;
    for (int i=0;i<like_trans.size();i++) {
      for (int j=0;j<like_trans[i].size();j++) {
        from[k] = pdat[ordered[i]].code;
        to[k] = pdat[ordered[j]].code;
		from_hcw[k] = pdat[ordered[i]].hcw;
		to_hcw[k] = pdat[ordered[j]].hcw;
		orderi[k] = ordered[i];
		orderj[k] = ordered[j];
        if (p.ali_file.compare("")==0) {
          likelihood[k] = like_trans[ordered[i]][ordered[j]].ns_lL_tot;
          consistency[k] = ThresholdsNSR(p, like_trans[ordered[i]][ordered[j]].ns_lL_tot);
          under_threshold[k] = (likelihood[k] < p.thresholdns);
        } else {
          likelihood[k] = like_trans[ordered[i]][ordered[j]].lL_tot;
          consistency[k] = ThresholdsR(p, like_trans[ordered[i]][ordered[j]].lL_tot);
          under_threshold[k] = (likelihood[k] < p.threshold);
        }
        ++k;
      }
    }
    return DataFrame::create( Named("from") = from,
                              Named("to")  = to,
                              Named("from_hcw")  = from_hcw,
                              Named("to_hcw")  = to_hcw,
                              Named("ordered_i")  = orderi,
                              Named("ordered_j")  = orderj,
                              Named("likelihood") = likelihood,
                              Named("consistency") = consistency,
                              Named("under_threshold") = under_threshold );
}

void LikelihoodOutput (run_params p, const vector<int>& ordered, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans) {
	for (int i=0;i<like_trans.size();i++) {
		for (int j=0;j<like_trans[i].size();j++) {
			Rcout << "From " << pdat[ordered[i]].code << " to " << pdat[ordered[j]].code << " ";
			if (p.ali_file.compare("")==0) {
				Rcout << like_trans[ordered[i]][ordered[j]].ns_lL_tot << " ";
				ThresholdsNS(p,like_trans[ordered[i]][ordered[j]].ns_lL_tot);
				Rcout << "\n";
			} else {
				Rcout << like_trans[ordered[i]][ordered[j]].lL_tot << " ";
				Thresholds(p,like_trans[ordered[i]][ordered[j]].lL_tot);
				Rcout << "\n";
			}
		}
	}
}

void FinalOutput (run_params p, const vector<string>& removed, const vector<string>& fixed) {
	if (removed.size()>0) {
		Rcout << "Individuals were removed due to missing or poor quality sequence data:\n";
		for (int i=0;i<removed.size();i++) {
			Rcout << removed[i] << "\n";
		}
	}
	if (fixed.size()>0) {
		Rcout << "Default location parameters were used to assess transmission events involving some individuals:\n";
		for (int i=0;i<fixed.size();i++) {
			Rcout << fixed[i] << "\n";
		}
	}
	Rcout << "Note that default parameters are biased towards individuals being present on the primary ward in question\n";
}
