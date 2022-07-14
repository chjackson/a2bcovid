using namespace std;
#include "basicmodel.h"

#include "Rcpp.h"

void FindVariants (vector<sparseseq>& variants, string& consensus, vector<pat>& pdat) {
    Rcpp::Function msg("message");
	msg("Find variants");
	for (int i=0;i<pdat.size();i++) {
        sparseseq s;
        for (int pos=0;pos<pdat[i].seq.size();pos++) {
			if (pdat[i].seq.compare(pos,1,consensus,pos,1)!=0) {
				if (pdat[i].seq.compare(pos,1,"A")==0||pdat[i].seq.compare(pos,1,"C")==0||pdat[i].seq.compare(pos,1,"G")==0||pdat[i].seq.compare(pos,1,"T")==0) {
					//Rcpp::Rcout << "Found variant " << pdat[i].code_match << " " << pos << " " << consensus[pos] << " " << pdat[i].seq[pos] << "\n";
					if (pos!=28827&&pos!=28828) {
						s.locus.push_back(pos);
						s.allele.push_back(pdat[i].seq[pos]);
					}
                }
            }
        }
        variants.push_back(s);
    }
}

void FindVariantsNoSeq (vector<sparseseq>& variants, vector<pat>& pdat) {
    Rcpp::Function msg("message");
	msg("Find variants NoSeq");
	for (int i=0;i<pdat.size();i++) {
        sparseseq s;
        variants.push_back(s);
    }
}


bool compare_allele(allele v1, allele v2) {
	return (v1.loc < v2.loc);
}

void ListAllVariantPositions (const vector<sparseseq>& variants, vector<allele>& allvar) {
	//Put into the vector allvar
	for (int i=0;i<variants.size();i++) {
		for (int j=0;j<variants[i].locus.size();j++) {
			allele a;
			a.loc=variants[i].locus[j];
			a.nuc=variants[i].allele[j];
			allvar.push_back(a);
		}
	}
<<<<<<< HEAD
    if (allvar.size()>=2) {
=======
    if (allvar.size()>1) {
>>>>>>> 64f475fd884a73fbd2c44440ccfdf53ea24e4771
        sort(allvar.begin(),allvar.end(),compare_allele);
        vector<int> to_rem;
        for (int i=1;i<allvar.size();i++) {
            if (allvar[i].loc==allvar[i-1].loc&&allvar[i].nuc==allvar[i-1].nuc) {
                to_rem.push_back(i);
            }
        }
        vector<int> keep;
        int index=0;
        for (int i=0;i<allvar.size();i++) {
            if (i==to_rem[index]) {
                index++;
            } else {
                keep.push_back(i);
            }
        }
        vector<allele> allvar_new;
        for (int i=0;i<keep.size();i++) {
            allvar_new.push_back(allvar[keep[i]]);
        }
        allvar=allvar_new;
    }
}


void FindAmbiguousVarPositions (const vector<allele>& allvar, vector<pat>& pdat, vector<int>& nloc_count) {
	//Variant positions for which at least one sequence has an N
	for (int i=0;i<pdat.size();i++) {
		int nc=0;
		for (int j=0;j<allvar.size();j++) {
			if (allvar[j].loc>pdat[i].seq.size()) {
				idat id;
				id.loc=allvar[j].loc;
				id.q=0;
				pdat[i].seq_uncertainty.push_back(id);
				nc++;
			} else {
				if (pdat[i].seq.compare(allvar[j].loc,1,"N")==0) {
					idat id;
					id.loc=allvar[j].loc;
					id.q=0;
					pdat[i].seq_uncertainty.push_back(id);
					nc++;
				}
			}
		}
		nloc_count.push_back(nc);
	}
}






