using namespace std;
#include "basicmodel.h"
#include "distributions.h"
#include "likelihoods.h"

double CalcLikelihood (run_params p, int sa, int sb, int da, int db, int eva, int evb, const vector<int> possible_tab) {
	double lL=0;
	double lL_seq=0;
	double lL_nonseq=0;
	//cout << sa << " " << sb << " " << da << " " << db << " " << eva << " " << evb << " ";
	for (int t=0;t<possible_tab.size();t++) {
		double L=0;

		if (db<possible_tab[t]) {  //Can't sequence B before transmission
			L=L-1e10;
		} else {
			//Sequence likelihood
			if (da<possible_tab[t]) {
				int evo_time=db-da;
				double evo=evo_time*p.rate;
				L=L+Poisson(eva,(p.seq_noise/2)); //D1 data - any novel variants must arise from noise
				L=L+Poisson(evb,evo+(p.seq_noise/2)); //Evolution between D1 and D2
			} else {
				int evo_time1=da-possible_tab[t];
				int evo_time2=db-possible_tab[t];
				double evo1=evo_time1*p.rate;
				double evo2=evo_time2*p.rate;
				L=L+Poisson(eva,evo1+(p.seq_noise/2)); //Evolution from transmission to D1
				L=L+Poisson(evb,evo2+(p.seq_noise/2)); //Evolution from transmission to D2
			}
		}
		lL_seq=lL_seq+exp(L);
		double lL_ns=0;
		double OGL=OffsetGammaCDFFlex(possible_tab[t]-sa,p.pa,p.pb,p.po); //Time from symptom onset of i to transmission
		L=L+OGL;
		lL_ns=lL_ns+OGL;
		double LN=LogNormal(sb-possible_tab[t],p.smu,p.ssigma);//Time from transmission to j becoming symptomatic
		L=L+LN;
		lL_ns=lL_ns+LN;
		lL_nonseq=lL_nonseq+exp(lL_ns);
		lL=lL+exp(L);
	}
	lL_seq=lL_seq/possible_tab.size();
	lL_seq=log(lL_seq);
	lL_nonseq=log(lL_nonseq);
	lL=log(lL);
	//cout << lL_seq << " " << lL_nonseq << " " << lL << "\n";
	return lL;
}
