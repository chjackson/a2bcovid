
#include "basicmodel.h"
#include "distributions.h"
#include "likelihoods.h"
#include "io.h"

using namespace std;
#include "Rcpp.h"

//Used for sorting data
bool comparetprob (tprob v1, tprob v2) {
	return (v1.time < v2.time);
}

void CalculateTDLikelihoods (run_params p, const vector<pat>& pdat, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, vector< vector<ijlike> >& like_trans) {
//	Rcpp::Rcout << "Calculate likelihoods\n";

Rcpp::Function msg("message");
  msg("Calculate likelihoods");
	for (int i=0;i<pdat.size();i++) {
		vector<ijlike> lt;
		for (int j=0;j<pdat.size();j++) {
			ijlike lij;
			lij.lL_tot=0;
			lij.ns_lL_tot=0;
			for (int k=0;k<pdat[i].time_s.time.size();k++) {
				lij.da.push_back(pdat[i].time_seq-pdat[i].time_s.time[k]);
			}
			for (int k=0;k<pdat[j].time_s.time.size();k++) {
				lij.db.push_back(pdat[j].time_seq-pdat[j].time_s.time[k]);
			}
			vector<tprob> contact_times_probs;
			if (pdat[i].code==pdat[j].code) { //Can't transmit to yourself
				lij.lL_tot=-1e10;
				lij.ns_lL_tot=-1e10;
				lij.min=0;
				lij.max=-1;
			} else {
				for (int k=0;k<pdat[i].locat.size();k++) {
					for (int l=0;l<pdat[j].locat.size();l++) {
						if (pdat[i].locat[k].ward.compare(pdat[j].locat[l].ward)==0&&pdat[i].locat[k].date==pdat[j].locat[l].date) {
							//Same time and place
							tprob t;
							t.time=pdat[i].locat[k].date;
							t.weight=pdat[i].locat[k].prob*pdat[j].locat[l].prob;
							contact_times_probs.push_back(t);
						}
					}
				}
				//Rcpp::Rcout << "Number of positive contacts " << contact_times_probs.size() << "\n";
				if (contact_times_probs.size()==0) {
					lij.lL_tot=-1e10;
					lij.ns_lL_tot=-1e10;
					lij.min=0;
					lij.max=-1;
				} else {
					RemoveDuplicateTimes(contact_times_probs);
					lij.min=contact_times_probs[0].time;
					lij.max=contact_times_probs[contact_times_probs.size()-1].time;
					//Rcpp::Rcout << lij.min << " " << lij.max << "\n";

					FillTimes(contact_times_probs);
					for (int k=0;k<contact_times_probs.size();k++) {
						double lL=LikelihoodFromItoJTimeK (i,j,k,p,contact_times_probs,seqdists,seqdists_c,pdat);
						lij.lL_tot=lij.lL_tot+exp(lL);
						lij.contact_times.push_back(contact_times_probs[k].time);
						lij.contact_likes.push_back(lL);
						double lLnoseq=NoSeqLikelihoodFromItoJTimeK (i,j,k,p,contact_times_probs,pdat);
						lij.noseq_likes.push_back(lLnoseq);
						lij.ns_lL_tot=lij.ns_lL_tot+exp(lLnoseq);
					}

					if (lij.lL_tot>0) {
						lij.lL_tot=log(lij.lL_tot);
					} else {
						lij.lL_tot=-1e10;
					}
					if (lij.ns_lL_tot>0) {
						lij.ns_lL_tot=log(lij.ns_lL_tot);
					} else {
						lij.ns_lL_tot=-1e10;
					}
				}
			}
			//Standardise impossible tranmissions - use this in optimising times to skip possibilities
			if (lij.lL_tot<-1e9) {
				lij.lL_tot=-1e10;
			}
			lt.push_back(lij);
		}
		like_trans.push_back(lt);
	}
}

void RemoveDuplicateTimes (vector<tprob>& contact_times_probs) {
	sort(contact_times_probs.begin(),contact_times_probs.end(),comparetprob);
	vector<int> to_rem;
	for (int k=0;k<contact_times_probs.size()-1;k++) {
		if (contact_times_probs[k].time==contact_times_probs[k+1].time) {
			if (contact_times_probs[k].weight<contact_times_probs[k+1].weight) {
				to_rem.push_back(k);
				contact_times_probs[k+1].weight = 1-((1-contact_times_probs[k+1].weight)*(1-contact_times_probs[k].weight));
			} else {
				to_rem.push_back(k+1);
				contact_times_probs[k].weight = 1-((1-contact_times_probs[k+1].weight)*(1-contact_times_probs[k].weight));
			}
		}
	}
	sort(to_rem.begin(),to_rem.end());
	reverse(to_rem.begin(),to_rem.end());
	for (int k=0;k<to_rem.size();k++) {
		contact_times_probs.erase(contact_times_probs.begin()+to_rem[k]);
	}
}

void FillTimes (vector<tprob>& contact_times_probs) {
	int size=contact_times_probs.size();
	int pos=0;
	while (pos<size-1) {
		if (contact_times_probs[pos+1].time!=contact_times_probs[pos].time+1) {
			tprob t;
			t.time=contact_times_probs[pos].time+1;
			t.weight=0;
			contact_times_probs.insert(contact_times_probs.begin()+pos+1,t);
			size++;
		}
		pos++;
	}
}

//New routine here.
double LikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector< vector<int> >& seqdists, const vector< vector<tpair> >& seqdists_c, const vector<pat>& pdat) {
	//Transmission from i to j at time contact_times[k].  Integrate sequence evolution date with key likelihood.
	double L=0;
	if (pdat[i].seq.size()==0||pdat[j].seq.size()==0) {
		L=-1e10;
	} else if (pdat[j].type==1) {//Flagged as community infection - can't be from i to j
		L=-1e10;
	} else {
		//Sequence component of the likelihoods
		if (pdat[j].time_seq<contact_times_probs[k].time) { //Data must be collected from j after the time of transmission
			L=L-1e10; //Can't transmit after collecting data from recipient
		} else {
			if (p.ali_file.compare("")!=0) {
				//Two sequence likelihoods -
				if (pdat[i].time_seq<contact_times_probs[k].time) {//First, case where transmission is after collecting the first sequence sample
					int evo_time=(pdat[j].time_seq-pdat[i].time_seq);//Amount of time over which evolution has occurred
					double evo=evo_time*p.rate;
					if (seqdists[i][j]==-1) {
						L=L+Poisson(9,evo+p.seq_noise);
					} else {
						L=L+Poisson(seqdists_c[i][j].from,(p.seq_noise/2)); //D1 data - any novel variants must arise from noise
						L=L+Poisson(seqdists_c[i][j].to,evo+(p.seq_noise/2)); //Evolution between D1 and D2
					}
				} else { //Second, case where transmission is before collecting the first sequence sample
					//Calculate times
					int evo_time1=pdat[i].time_seq-contact_times_probs[k].time;
					int evo_time2=pdat[j].time_seq-contact_times_probs[k].time;
					double evo1=evo_time1*p.rate;
					double evo2=evo_time2*p.rate;
					if (seqdists[i][j]==-1) {
						L=L+Poisson(9,evo1+(p.seq_noise/2));
						L=L+Poisson(9,evo2+(p.seq_noise/2));
					} else {
						L=L+Poisson(seqdists_c[i][j].from,evo1+(p.seq_noise/2)); //Evolution from transmission to D1
						L=L+Poisson(seqdists_c[i][j].to,evo2+(p.seq_noise/2)); //Evolution from transmission to D2
					}
				}
			}
			double Ladd=0;
			for (int q=0;q<pdat[i].time_s.prob.size();q++) {
				double OGL=OffsetGammaCDFFlex(contact_times_probs[k].time-pdat[i].time_s.time[q],p.pa,p.pb,p.po); //Time from symptom onset of i to transmission
				Ladd=Ladd+(pdat[i].time_s.prob[q]*exp(OGL));
				if (pdat[i].time_s.time.size()>1) {
				}

			}
			Ladd=log(Ladd);
			L=L+Ladd;
			if (contact_times_probs[k].weight==0) {
				L=L-1e10; //Can't transmit when there is no contact between individuals
			} else {
				L=L+log(contact_times_probs[k].weight); //Probability of contact between individuals
			}
			Ladd=0;
			for (int q=0;q<pdat[i].time_s.prob.size();q++) {
				double LN=LogNormal(pdat[j].time_s.time[q]-contact_times_probs[k].time,p.smu,p.ssigma);//Time from transmission to j becoming symptomatic
				Ladd=Ladd+(pdat[i].time_s.prob[q]*exp(LN));
			}
			Ladd=log(Ladd);
			L=L+Ladd;
		}
	}
	return L;
}


double NoSeqLikelihoodFromItoJTimeK (int i, int j, int k, const run_params& p, const vector<tprob>& contact_times_probs, const vector<pat>& pdat) {
	//Transmission from i to j at time contact_times[k].  Integrate sequence evolution date with key likelihood.
	double L=0;
	double Ladd=0;
	for (int q=0;q<pdat[i].time_s.prob.size();q++) {
		double OGL=OffsetGammaCDFFlex(contact_times_probs[k].time-pdat[i].time_s.time[q],p.pa,p.pb,p.po); //Time from symptom onset of i to transmission
		Ladd=Ladd+(pdat[i].time_s.prob[q]*exp(OGL));
	}
	Ladd=log(Ladd);
	L=L+Ladd;
	if (contact_times_probs[k].weight==0) {
		L=L-1e10; //Can't transmit when there is no contact between individuals
	} else {
		L=L+log(contact_times_probs[k].weight); //Probability of contact between individuals
	}
	Ladd=0;
	for (int q=0;q<pdat[i].time_s.prob.size();q++) {
		double LN=LogNormal(pdat[j].time_s.time[q]-contact_times_probs[k].time,p.smu,p.ssigma);//Time from transmission to j becoming symptomatic
		Ladd=Ladd+(pdat[i].time_s.prob[q]*exp(LN));
	}
	Ladd=log(Ladd);
	L=L+Ladd;
	//Rcpp::Rcout << "Log " << L << "\n";
	return L;
}

/*
//This code isn't currently used, but could be of use if we want to divide the likelihoods into subsets later.
void FindLikelihoodSubsetsIJ(run_params p, const vector< vector<ijlike> >& like_trans, vector< vector<int> >& subsets) {
	vector <int> range;
	for (int i=0;i<like_trans.size();i++) {
		range.push_back(i);
	}
	while (range.size()>0) {
		vector<int> sset;
		sset.clear();
		sset.push_back(range[0]);
		int add=1;
		while (add==1) {
			int ss=sset.size();
			add=0;
			for (int i=0;i<ss;i++) {
				for (int j=0;j<like_trans.size();j++) {
					if (like_trans[sset[i]][j].lL_tot>p.threshold) {
						//Rcpp::Rcout << "Add " << sset[i] << " " << j << "\n";
						sset.push_back(j);
					}
					if (like_trans[j][sset[i]].lL_tot>p.threshold) {
						//Rcpp::Rcout << "Add " << sset[i] << " " << j << "\n";
						sset.push_back(j);
					}
				}
			}
			sort(sset.begin(),sset.end());
			sset.erase(unique(sset.begin(),sset.end()),sset.end());
			//Rcpp::Rcout << sset.size() << "\n";
			if (sset.size()>ss) {
				add=1;
			}
		}
		//Add sset to subsets
		//Rcpp::Rcout << "Push back to subsets\n";
		subsets.push_back(sset);
		for (int i=0;i<sset.size();i++) {
			for (int j=0;j<range.size();j++) {
				if (range[j]==sset[i]) {
					range.erase(range.begin()+j);
					break;
				}
			}
		}
	}
}
*/


Rcpp::String ThresholdsR (run_params p, int i, int j, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans) {
  Rcpp::String out;
	double t95=0;
	double t99=0;
	int error=0;
	for (int k1=0;k1<pdat[i].time_s.time.size();k1++) {
		if (error==0) {
			int da=like_trans[i][j].da[k1]+10;
			for (int k2=0;k2<pdat[j].time_s.time.size();k2++) {
				if (error==0) {
					int db=like_trans[i][j].db[k2]+10;
					if (da<0||db<0||da>51||db>51) {
						out = "Unlikely";
						error=1;
						break;
					}
					t95=t95+pdat[i].time_s.prob[k1]*pdat[j].time_s.prob[k2]*p.threshold95[da][db];
					t99=t99+pdat[i].time_s.prob[k1]*pdat[j].time_s.prob[k2]*p.threshold99[da][db];
				}
			}
		}
	}
	if (error==0) {
		if (like_trans[i][j].lL_tot>t95) {
			out = "Consistent"; //p>0.05
		} else if (like_trans[i][j].lL_tot>t99) {
			out = "Borderline"; //0.05>p>0.01
		} else {
			out = "Unlikely"; //0.01>p
		}
	}
  return out;
}

Rcpp::String ThresholdsNSR (run_params p, double L) {
	Rcpp::String out;
	if (L>p.t95NS) {
		out = "Consistent"; //p>0.05
	} else if (L>p.t99NS) {
		out = "Borderline"; //0.05>p>0.01
	} else {
		out = "Unlikely";
	}
	return out;
}

void Thresholds (run_params p, int i, int j, const vector<pat>& pdat, const vector< vector<ijlike> >& like_trans) {
	//Calculate threshold; Requires probabilistic sum over D_i - S_i and D_j - S_j
	double t95=0;
	double t99=0;
	int error=0;
	for (int k1=0;k1<pdat[i].time_s.time.size();k1++) {
		if (error==0) {
			int da=like_trans[i][j].da[k1]+10;
			for (int k2=0;k2<pdat[j].time_s.time.size();k2++) {
				if (error==0) {
					int db=like_trans[i][j].db[k2]+10;
					if (da<0||db<0||da>51||db>51) {
						Rcpp::Rcout << "Unlikely ";
						error=1;
						break;
					}
					t95=t95+pdat[i].time_s.prob[k1]*pdat[j].time_s.prob[k2]*p.threshold95[da][db];
					t99=t99+pdat[i].time_s.prob[k1]*pdat[j].time_s.prob[k2]*p.threshold99[da][db];
				}
			}
		}
	}
	if (error==0) {
		if (like_trans[i][j].lL_tot>t95) {
			Rcpp::Rcout << "Consistent "; //p>0.05
		} else if (like_trans[i][j].lL_tot>t99) {
			Rcpp::Rcout << "Borderline "; //0.05>p>0.01
		} else {
			Rcpp::Rcout << "Unlikely "; //0.01>p
		}
	}
}

void ThresholdsNS (run_params p, double L) {
	if (L>p.t95NS) {
		Rcpp::Rcout << "Consistent "; //p>0.05
	} else if (L>p.t99NS) {
		Rcpp::Rcout << "Borderline "; //0.05>p>0.01
	} else {
		Rcpp::Rcout << "Unlikely "; //0.01>p
	}
}
