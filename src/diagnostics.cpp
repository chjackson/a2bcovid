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

DataFrame LikelihoodOutputR(run_params p, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans) {
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
    NumericVector likelihood(n);
    CharacterVector consistency(n);
    LogicalVector under_threshold(n);
    k = 0;
    for (int i=0;i<like_trans.size();i++) {
      for (int j=0;j<like_trans[i].size();j++) {
        from[k] = pdat[i].code;
        to[k] = pdat[j].code;
        if (p.data_type==0) {
          likelihood[k] = like_trans[i][j].ns_lL_tot;
          consistency[k] = ThresholdsNSR(p, like_trans[i][j].ns_lL_tot);
          under_threshold[k] = (likelihood[k] < p.thresholdns);
        } else {
          likelihood[k] = like_trans[i][j].lL_tot;
          consistency[k] = ThresholdsR(p, like_trans[i][j].lL_tot);
          under_threshold[k] = (likelihood[k] < p.threshold);
        }
        ++k;
      }
    }
    return DataFrame::create( Named("from") = from,
                              Named("to")  = to,
                              Named("likelihood") = likelihood,
                              Named("consistency") = consistency,
                              Named("under_threshold") = under_threshold );
}

void LikelihoodOutput (run_params p, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans) {
	for (int i=0;i<like_trans.size();i++) {
		for (int j=0;j<like_trans[i].size();j++) {
			Rcout << "From " << pdat[i].code << " to " << pdat[j].code << " ";
			if (p.data_type==0) {
				Rcout << like_trans[i][j].ns_lL_tot << " ";
				ThresholdsNS(p,like_trans[i][j].ns_lL_tot);
				Rcout << "\n";
			} else {
				Rcout << like_trans[i][j].lL_tot << " ";
				Thresholds(p,like_trans[i][j].lL_tot);
				Rcout << "\n";
			}
		}
	}
	Rcout << p.data_type << "\n";
}

void FinalOutput (run_params p, const vector<string>& removed, const vector<string>& fixed) {
	if (removed.size()>0) {
		Rcout << "Individuals were removed due to missing or poor quality sequence data:\n";
		for (int i=0;i<removed.size();i++) {
			Rcout << removed[i] << "\n";
		}
	}
	Rcout << "Can rerun with --data_type 0 to perform an assessment based on symptom and location data alone\n\n";

	if (fixed.size()>0) {
		Rcout << "Default location parameters were used to assess transmission events involving some individuals:\n";
		for (int i=0;i<fixed.size();i++) {
			Rcout << fixed[i] << "\n";
		}
	}
	Rcout << "Note that default parameters are biased towards individuals being present on the primary ward in question\n";


	if (p.data_type==1) {
		Rcout << "Can rerun with --data_type 2 or --data_type 3 if location data are available\n\n";
	}
	if (p.data_type==2) {
		Rcout << "Can rerun with --data_type 3 if HCW location data are available\n\n";
	}

}
