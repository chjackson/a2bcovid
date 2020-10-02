#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "basicmodel.h"
#include "utilities.h"
#include "distributions.h"
#include "likelihoods.h"
#include "io.h"

#include "Rcpp.h"

void RemoveDuplicatesNoSeq (vector<pat>& pdat) {
	vector<int> to_rem;
	for (int i=1;i<pdat.size();i++) {
		if (pdat[i].code==pdat[i-1].code) {
			to_rem.push_back(i);
		}
	}
	sort(to_rem.begin(),to_rem.end());
	reverse(to_rem.begin(),to_rem.end());
	for (int i=0;i<to_rem.size();i++) {
		pdat.erase(pdat.begin()+to_rem[i]);
	}
}

void CorrectNames (vector<string>& names) {
	for (int i=0;i<names.size();i++) {
		names[i].erase(0,1);
		stringstream a(names[i]);
		string segment;
		vector<string> seglist;
		while(getline(a, segment, '/')) {
		   seglist.push_back(segment);
		}
		names[i]=seglist[0];
	}
}



void MatchSequencetoPatient (vector<pat>& pdat) {
	ifstream match_file;
	match_file.open("../Model_data/Sequence_matches.dat");
	for (int i=0;i<1000000;i++) {
		string code;
		string code_match;
		if (!(match_file >> code)) break;
		if (!(match_file >> code_match)) break;
		for (int j=0;j<pdat.size();j++) {
			//Rcpp::Rcout << pdat[i].code << " " << code << "\n";
			if (pdat[j].code==code) {
			//	Rcpp::Rcout << "Match\n";
				pdat[j].code_match=code_match;
			}
		}
	}
}

void IncorporateSequenceData (run_params p, vector<pat>& pdat, vector<string>& names, vector<string>& seqs) {
	for (int i=0;i<names.size();i++) {
		//Rcpp::Rcout << "Name " << i << " of " << names.size() << " " << names[i] << "\n";
		for (int j=0;j<pdat.size();j++) {
			if (pdat[j].code_match.compare(names[i])==0) {
				if (p.diagnostic==1) {
					Rcpp::Rcout << "Match " << j << " " << pdat[j].code << " " << pdat[j].code_match << "\n";
				}
				pdat[j].seq=seqs[i];
			}
		}
	}
}

void CountNs (vector<pat>& pdat) {
	//Number of Ns in each sequence.  Tiebreaker for samples when there are more than one from the same day.
	for (int i=0;i<pdat.size();i++) {
		int nn=0;
		//cout << i << " " << pdat[i].seq.size() << "\n";
		for (int j=0;j<pdat[i].seq.size();j++) {
			if (pdat[i].seq.compare(j,1,"N")==0) {
				nn++;
			}
		}
		//cout << i << " " << pdat[i].seq.size() << " " << nn << "\n";
		pdat[i].seq_n=nn;
	}
}

void FindConsensus (string& consensus, vector<string>& seqs) {
	consensus=seqs[0];
	Rcpp::Rcout << "Find consensus of all input sequences\n";
    int nA=0;
    int nC=0;
    int nG=0;
    int nT=0;
    for (int pos=0;pos<seqs[0].size();pos++) {
        nA=0;
        nC=0;
        nG=0;
        nT=0;
        for (int seq=0;seq<seqs.size();seq++) {
            if (seqs[seq][pos]=='A') {
                nA++;
            }
            if (seqs[seq][pos]=='C') {
                nC++;
            }
            if (seqs[seq][pos]=='G') {
                nG++;
            }
            if (seqs[seq][pos]=='T') {
                nT++;
            }
        }
        int max=nA;
        consensus[pos]='A';
        if (nC>max) {
            max=nC;
            consensus[pos]='C';
        }
        if (nG>max) {
            max=nG;
            consensus[pos]='G';
        }
        if (nT>max) {
            consensus[pos]='T';
            max=nT;
        }
        if (max==0) {
            consensus[pos]='-';
        }
    }
}

void FindSequenceLength(int& seq_len, string& consensus, vector<pat> pdat) {
	for (int i=0;i<pdat.size();i++) {
		if (seq_len<pdat[i].seq.size()) {
			seq_len=pdat[i].seq.size();
			consensus=pdat[i].seq;
		}
	}
}

void FindPairwiseDistances (run_params p, vector< vector<int> >& seqdists, vector<sparseseq>& variants, vector<pat>& pdat) {
	vector<int> zeros(pdat.size(),0);
	for (int i=0;i<pdat.size();i++) {
		seqdists.push_back(zeros);
	}
    for (int i=0;i<pdat.size();i++) {
        for (int j=i+1;j<pdat.size();j++) {
            int dist=0;
            //Find unique difference positions;
            vector<int> uniq;
            for (int k=0;k<variants[i].locus.size();k++) {
                uniq.push_back(variants[i].locus[k]);
            }
            for (int k=0;k<variants[j].locus.size();k++) {
                uniq.push_back(variants[j].locus[k]);
            }
            sort(uniq.begin(),uniq.end());
            uniq.erase(unique(uniq.begin(),uniq.end()),uniq.end());
            for (int k=0;k<uniq.size();k++) {
                if (pdat[i].seq[uniq[k]]=='A'||pdat[i].seq[uniq[k]]=='C'||pdat[i].seq[uniq[k]]=='G'||pdat[i].seq[uniq[k]]=='T') {
                    if (pdat[j].seq[uniq[k]]=='A'||pdat[j].seq[uniq[k]]=='C'||pdat[j].seq[uniq[k]]=='G'||pdat[j].seq[uniq[k]]=='T') {
                        if (pdat[i].seq[uniq[k]]!=pdat[j].seq[uniq[k]]) {
                            dist++;
                        }
                    }
                }
            }
			seqdists[i][j]=dist;
			seqdists[j][i]=dist;
        }
    }
	for (int i=0;i<pdat.size();i++) {
		if (pdat[i].seq.size()==0) {
			for (int j=0;j<pdat.size();j++){
				seqdists[i][j]=-1;
				seqdists[j][i]=-1;
			}
			seqdists[i][i]=0;
		}
	}
}

void FindPairwiseDistancesNoSeq (run_params p, vector< vector<int> >& seqdists, vector<sparseseq>& variants, vector<pat>& pdat) {
	vector<int> zeros(pdat.size(),0);
	for (int i=0;i<pdat.size();i++) {
		seqdists.push_back(zeros);
	}
    for (int i=0;i<pdat.size();i++) {
        for (int j=i+1;j<pdat.size();j++) {
           	seqdists[i][j]=0;
			seqdists[j][i]=0;
        }
    }
}


void FromConsensusDistances (const vector<sparseseq>& variants, vector< vector<tpair> >& seqdists_c) {
	//Construct a pairwise consensus between i and j as having a list of variants present in both sequences.
	//Find the number of variants each sequence would need to gain in order to be constructed from that consensus.
	for (int i=0;i<variants.size();i++) {
		vector<tpair> vc;
		for (int j=0;j<variants.size();j++) {
			tpair p;
			p.from=0;
			p.to=0;
			//We want the number of variants in i not in j and vice versa.
			for (int k=0;k<variants[i].locus.size();k++) {
				int mycount = count(variants[j].locus.begin(), variants[j].locus.end(), variants[i].locus[k]);
				if (mycount==0) {
					p.from++;
				}
			}
			for (int k=0;k<variants[j].locus.size();k++) {
				int mycount = count(variants[i].locus.begin(), variants[i].locus.end(), variants[j].locus[k]);
				if (mycount==0) {
					p.to++;
				}
			}
			vc.push_back(p);
		}
		seqdists_c.push_back(vc);
	}
}

void FromConsensusDistancesNoSeq (const vector<sparseseq>& variants, vector< vector<tpair> >& seqdists_c) {
	//Construct a pairwise consensus between i and j as having a list of variants present in both sequences.
	//Find the number of variants each sequence would need to gain in order to be constructed from that consensus.
	for (int i=0;i<variants.size();i++) {
		vector<tpair> vc;
		for (int j=0;j<variants.size();j++) {
			tpair p;
			p.from=0;
			p.to=0;
			vc.push_back(p);
		}
		seqdists_c.push_back(vc);
	}
}


void ProcessSymptomUncertainty (run_params& p, vector<pat>& pdat) {
	//Deals with cases where we have positive test dates but not dates of becoming symptomatic
	//Find samples with uncertain symptomatic times
	vector<int> uncertain_times;
	for (int i=0;i<pdat.size();i++) {
		if (pdat[i].time_s_cert==0) {
			uncertain_times.push_back(i);
		}
	}
	if (uncertain_times.size()>0) {
		Rcpp::Rcout << "Estimate unknown times of symptom onset using CUH mean statistics\n";
		for (int i=0;i<pdat.size();i++) {
			if (pdat[i].time_s_cert==0) {
				pdat[i].time_s=pdat[i].time_s-floor(p.uct_mean+0.5);
			}
		}
	}
}

void RemoveIndividualsNoSequence(const run_params p, vector<pat>& pdat, vector<string>& removed) {
	vector<int> to_rem;
	for (int i=0;i<pdat.size();i++) {
		if (pdat[i].seq.size()==0) {
			if (p.diagnostic==1) {
				Rcpp::Rcout << "No sequence data for individual " << pdat[i].code << " " << pdat[i].code_match << ": Excluding from input\n";
			}
			to_rem.push_back(i);
		}
	}
	sort(to_rem.begin(),to_rem.end());
	reverse(to_rem.begin(),to_rem.end());
	for (int i=0;i<to_rem.size();i++) {
		removed.push_back(pdat[to_rem[i]].code);
		pdat.erase(pdat.begin()+to_rem[i]);
	}
}

void RemoveIndividualsSequenceQuality(const run_params p, vector<pat>& pdat, vector<string>& removed) {
	vector<int> to_rem;
	double maxe=1-p.min_qual;
	for (int i=0;i<pdat.size();i++) {
		if (pdat[i].seq_n>maxe*pdat[i].seq.size()) {
			if (p.diagnostic==1) {
				Rcpp::Rcout << "Sequence data for individual " << pdat[i].code << " has low quality (% of on-ambiguous sites): Excluding from input\n";
			}
			to_rem.push_back(i);
		}
	}
	if (to_rem.size()>0) {
		if (p.diagnostic==1) {
			Rcpp::Rcout << "Required quality can be set with the flag --min_qual\n";
		}
	}
	sort(to_rem.begin(),to_rem.end());
	reverse(to_rem.begin(),to_rem.end());
	for (int i=0;i<to_rem.size();i++) {
		removed.push_back(pdat[to_rem[i]].code);
		pdat.erase(pdat.begin()+to_rem[i]);
	}
}


void RemoveIndividualsMultipleN (const run_params p, vector<int>& nloc_count, vector<pat>& pdat, vector<sparseseq>& variants, vector<string>& removed) {
	vector<int> to_rem;
	for (int i=pdat.size()-1;i>=0;i--) {
		if (nloc_count[i]>p.max_n) {
			if (p.diagnostic==1) {
				Rcpp::Rcout << "Sequence data for individual " << pdat[i].code << " has too many ambiguous sites at variant loci: Excluding from input\n";
			}
			to_rem.push_back(i);
		}
	}
	if (to_rem.size()>0) {
		if (p.diagnostic==1) {
			Rcpp::Rcout << "Maximum number can be set with the flag --max_n\n";
		}
	}
	for (int i=0;i<to_rem.size();i++) {
		removed.push_back(pdat[i].code);
		pdat.erase(pdat.begin()+to_rem[i]);
		variants.erase(variants.begin()+to_rem[i]);
		nloc_count.erase(nloc_count.begin()+to_rem[i]);
	}
}

string FindMostCommonWard(vector<pat>& pdat) {
	vector<string> ward;
	for (int i=0;i<pdat.size();i++) {
		for (int j=0;j<pdat[i].locat.size();j++) {
			ward.push_back(pdat[i].locat[j].ward);
		}
	}
	string mostc;
	if (ward.size()>0) {
		vector<string> ward_u=ward;
		ward_u.erase(unique(ward_u.begin(),ward_u.end()),ward_u.end());
		vector<int> count;
		for (int i=0;i<ward_u.size();i++) {
			int c=std::count(ward.begin(),ward.end(),ward_u[i]);
			count.push_back(c);
		}
		int max=-1;
		for (int i=0;i<ward_u.size();i++) {
			if (count[i]>max) {
				mostc=ward_u[i];
				max=count[i];
			}
		}
	} else {
		mostc="WARD_01";
	}
	return mostc;
}

void FixIndivudualsNoLocation (const run_params p, vector<pat>& pdat, vector<string>& fixed) {
	//People are there by default according to their status
	string mostc=FindMostCommonWard(pdat);
//	Rcpp::Rcout << "Most common " << mostc << "\n";
	//Find minimum and maximum times
	int min=100000;
	int max=-100000;
	for (int i=0;i<pdat.size();i++) {
		for (int j=0;j<pdat[i].locat.size();j++) {
			if (pdat[i].locat[j].date>max) {
				max=pdat[i].locat[j].date;
			}
			if (pdat[i].locat[j].date<min) {
				min=pdat[i].locat[j].date;
			}
		}
		if (pdat[i].time_s+25>max) {
			max=pdat[i].time_s+25;
		}
		if (pdat[i].time_s-15<min) {
			min=pdat[i].time_s-15;
		}
	}
	//Rcpp::Rcout << max << " " << min << " " << p.hcw_location_default << " " << p.pat_location_default << "\n";
	//If no location data, assign probabilities by individual status
	for (int i=0;i<pdat.size();i++) {
		if (pdat[i].locat.size()==0) {
			fixed.push_back(pdat[i].code);
			loc l;
			l.ward=mostc;
			if (pdat[i].hcw==1) {
				l.prob=p.hcw_location_default;
			} else {
				l.prob=p.pat_location_default;
			}
			for (int j=min;j<max;j++) {
				l.date=j;
				pdat[i].locat.push_back(l);
			}
		}
	}
}

void RemoveRepeatPatients (run_params p, vector<pat>& pdat) {
	//Removes individuals with more than one listing i.e. sequence.  Takes the first sequence ech time.  If sequences have identical dates, takes the one with the fewest N nucleotides.  Use removed sequences to repair any gaps in the
	vector<int> rem;
	double maxe=1-p.min_qual;
	for (int i=0;i<pdat.size();i++) {
		if (pdat[i].seq_n<maxe*pdat[i].seq.size()) {
			//Remove anything later
			for (int j=0;j<pdat.size();j++) {
				if (pdat[i].code==pdat[j].code) {
					if (pdat[i].time_seq<pdat[j].time_seq) {
						rem.push_back(j);
						//Rcpp::Rcout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
						RepairSequence(i,j,pdat);
					}
					//Remove anything with the same time of worse quality
					if (pdat[i].time_seq==pdat[j].time_seq) {
						if (pdat[i].seq_n<pdat[j].seq_n) {
							rem.push_back(j);
							//Rcpp::Rcout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
							RepairSequence(i,j,pdat);
						}
					}
					//Remove anything of poor quality
					if (pdat[j].seq_n>maxe*pdat[i].seq.size()) {
						rem.push_back(j);
						//Rcpp::Rcout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
						RepairSequence(i,j,pdat);
					}
				}
			}
		} else {
			//Low quality.  Remove anything of worse quality
			for (int j=0;j<pdat.size();j++) {
				if (pdat[i].code==pdat[j].code) {
					if (pdat[i].seq_n<pdat[j].seq_n) {
						rem.push_back(j);
						//Rcpp::Rcout << "Remove " << pdat[j].code << " " << pdat[j].code_match << " " << j << " " << i << "\n";
						RepairSequence(i,j,pdat);
					}
				}
			}
		}
	}
	sort(rem.begin(),rem.end());
	rem.erase(unique(rem.begin(),rem.end()),rem.end());
	reverse(rem.begin(),rem.end());
	for (int i=0;i<rem.size();i++) {
		pdat.erase(pdat.begin()+rem[i]);
	}
	if (p.diagnostic==1) {
	Rcpp::Rcout << "Data quality\n";
		for (int i=0;i<pdat.size();i++) {
			Rcpp::Rcout << pdat[i].code << " " << pdat[i].code_match << " #N = " << pdat[i].seq_n << "\n";
		}
	}
}

void RepairSequence(int i, int j, vector<pat>& pdat) {
	for (int k=0;k<pdat[i].seq.size();k++) {
		if (pdat[i].seq.compare(k,1,"N")==0) {
			if (pdat[j].seq.compare(k,1,"N")!=0) {
				pdat[i].seq[k]=pdat[j].seq[k];
			}
		}
	}
}

void FindOrdering (const vector< vector<ijlike> >& like_trans, vector<int>& ordered) {
	vector< vector<int> > categories;
	FindLikelihoodCategories(like_trans,categories);
	vector< vector<int> > dists;
	FindCategoryDistances(categories,dists);
	vector< vector<int> > subsets;
	FindSubsets(dists,subsets);
	for (int i=0;i<subsets.size();i++) {
		for (int j=0;j<subsets[i].size();j++) {
			ordered.push_back(subsets[i][j]);
		}
	}
}

void FindLikelihoodCategories(const vector< vector<ijlike> >& like_trans, vector< vector<int> >& categories) {
	for (int i=0;i<like_trans.size();i++) {
		vector<int> c;
		for (int j=0;j<like_trans[i].size();j++) {
			if (like_trans[i][j].lL_tot>-8.15176) {
				c.push_back(1);
			} else if (like_trans[i][j].lL_tot>-10.1578) {
				c.push_back(2);
			} else {
				c.push_back(3);
			}
		}
		categories.push_back(c);
	}
}

void FindCategoryDistances(const vector< vector<int> >& categories, vector< vector<int> >& dists) {
	for (int i=0;i<categories.size();i++) {
		vector<int> c;
		for (int j=0;j<categories.size();j++) {
			int d=0;
			if (i!=j) {
				for (int k=0;k<categories[i].size();k++) {
					if (k!=i&&k!=j) {
						d=d+abs(categories[i][k]-categories[j][k]);
					}
					if (k==i) {
						d=d+categories[i][j];
					}
					if (k==j){
						d=d+categories[j][i];
					}
				}
			}
			c.push_back(d);
		}
		dists.push_back(c);
	}
}

void FindSubsets (const vector< vector<int> >& dists, vector< vector<int> >& subsets) {
	int done=0;
	while (done==0) {
		int min=100;
		int mi=-1;
		int mj=-1;
		for (int i=0;i<dists.size();i++) {
			for (int j=0;j<dists[i].size();j++) {
				if (i!=j&&dists[i][j]<min) {
					int found=0;
					for (int k=0;k<subsets.size();k++) {
						for (int l=0;l<subsets[k].size();l++) {
							if (subsets[k][l]==i) {
								found++;
							}
							if (subsets[k][l]==j) {
								found++;
							}
						}
					}
					if (found<2) {
						min=dists[i][j];
						mi=i;
						mj=j;
					}
				}
			}
		}
		//cout << mi << " " << mj << "\n";
		if (mi>=0&&mj>=0) {
			int placed=0;
			for (int k=0;k<subsets.size();k++) {
				for (int l=0;l<subsets[k].size();l++) {
					if (subsets[k][l]==mi) {
						subsets[k].push_back(mj);
						placed=1;
						break;
					}
					if (subsets[k][l]==mj) {
						subsets[k].push_back(mi);
						placed=1;
						break;
					}
				}
			}
			if (placed==0) {
				vector<int> n;
				n.push_back(mi);
				n.push_back(mj);
				subsets.push_back(n);
			}
		} else {
			done=1;
		}
	}
}
