
#include "basicmodel.h"
#include "thresholds.h"
#include "distributions.h"

#include "Rcpp.h"

using namespace std;

void CalculateThresholdsFull (run_params p) {
	Rcpp::Rcout << "Calculation for complete likelihood\n";
	vector<double> allstats;
	//Assume that S1 is at time zero
	double n=0.772469;
	double r=(0.0008*29782)/365.25;
	for (int s2=-11;s2<=87;s2++) {
		Rcpp::Rcout << "Track progress : S2=" << s2 << "\n";
		for (int d2=s2-5;d2<=s2+20;d2++) {
			for (int d1=-5;d1<=20;d1++){
				for (int h=0;h<=10;h++) {
					for (int b=0;b<=h;b++) {
						int h1=b;
						int h2=h-b;
						double logT=0;
						for (int t=-11;t<=16;t++) {
							if (t<=s2&&d2>=t) {
								if (d1<t) {
									int evotime2=d2-d1;
									double L=OffsetGammaCDFFlex(t,p.pa,p.pb,p.po);
									L=L+LogNormal(s2-t,p.smu,p.ssigma);
									double lambda=n/2;
									L=L+Poisson(h1,lambda);
									lambda=(n/2)+(evotime2*r);
									L=L+Poisson(h2,lambda);
									L=exp(L);
									logT=logT+L;
								} else {
									int evotime1=d1-t;
									int evotime2=d2-t;
									double L=OffsetGammaCDFFlex(t,p.pa,p.pb,p.po);
									L=L+LogNormal(s2-t,p.smu,p.ssigma);
									double lambda=(n/2)+(evotime1*r);
									L=L+Poisson(h1,lambda);
									lambda=(n/2)+(evotime2*r);
									L=L+Poisson(h2,lambda);
									L=exp(L);
									logT=logT+L;
								}
							}
						}
						logT=log(logT);
						allstats.push_back(logT);
					}
				}
			}
		}
	}
	Rcpp::Rcout << "Size " << allstats.size() << "\n";
	ofstream tf_file;
	tf_file.open("Threshold_data.out");
	for (int i=0;i<allstats.size();i++) {
		tf_file << allstats[i] << "\n";
	}
	//Calculation
	double tot=0;
	Rcpp::Rcout << "Do sum\n";
	for (int i=0;i<allstats.size();i++) {
		allstats[i]=exp(allstats[i]);
		tot=tot+allstats[i];
	}
	Rcpp::Rcout << "Do sort\n";
	sort(allstats.begin(),allstats.end());
	reverse(allstats.begin(),allstats.end());
	int found95=0;
	int found99=0;
	int found999=0;
	int found59=0;
	Rcpp::Rcout << "Do search\n";
	for (int i=1;i<allstats.size();i++) {
		double unlog=log(allstats[i]);
		allstats[i]=allstats[i]+allstats[i-1];
		if (allstats[i]>tot*0.95&&found95==0) {
			found95=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
		}
		if (allstats[i]>tot*0.99&&found99==0) {
			found99=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
		}
		if (allstats[i]>tot*0.999&&found999==0) {
			found999=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
		}
		if (allstats[i]>tot*0.99999&&found59==0) {
			found59=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
			break;
		}
	}
}


void CalculateThresholdsNoSeq (run_params p) {
	Rcpp::Rcout << "Calculation for likelihood with sequence information removed\n";
	vector<double> allstats;
	//Assume that S1 is at time zero
	stat ss;
	ss.s1=0;
	ss.h1=0;
	ss.h2=0;
	ss.h=0;
	ss.b=0;
	for (int s2=-11;s2<=87;s2++) {
		ss.s2=s2;
		double L=0;
		for (int t=-11;t<=16;t++) {
			if (t<=s2) {
				double Ln=OffsetGammaCDFFlex(t,p.pa,p.pb,p.po);
				Ln=Ln+LogNormal(s2-t,p.smu,p.ssigma);
				Ln=exp(Ln);
				L=L+Ln;
			}
		}
		L=log(L);
		ss.L=L;
		Rcpp::Rcout << s2 << " " << L << "\n";
		allstats.push_back(L);
	}
	Rcpp::Rcout << "Size " << allstats.size() << "\n";
	ofstream tf_file;
	tf_file.open("Threshold_data_noseq.out");
	for (int i=0;i<allstats.size();i++) {
		tf_file << allstats[i] << "\n";
	}
	//Calculation
	double tot=0;
	Rcpp::Rcout << "Do sum\n";
	for (int i=0;i<allstats.size();i++) {
		allstats[i]=exp(allstats[i]);
		tot=tot+allstats[i];
	}
	Rcpp::Rcout << "Do sort\n";
	sort(allstats.begin(),allstats.end());
	reverse(allstats.begin(),allstats.end());
	int found95=0;
	int found99=0;
	int found999=0;
	int found59=0;
	Rcpp::Rcout << "Do search\n";
	for (int i=1;i<allstats.size();i++) {
		double unlog=log(allstats[i]);
		allstats[i]=allstats[i]+allstats[i-1];
		if (allstats[i]>tot*0.95&&found95==0) {
			found95=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
		}
		if (allstats[i]>tot*0.99&&found99==0) {
			found99=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
		}
		if (allstats[i]>tot*0.999&&found999==0) {
			found999=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
		}
		if (allstats[i]>tot*0.99999&&found59==0) {
			found59=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << allstats[i]/tot << "\n";
			break;
		}
	}
}


void ReCalculateThresholds (run_params p) {
	ifstream tf_file;
	tf_file.open("Threshold_data.out");
	vector<double> L;
	double x;
	Rcpp::Rcout << "Reading data\n";
	for (int i=0;i<1000000000;i++) {
		if (!(tf_file >> x)) break;
		L.push_back(x);
	}
	Rcpp::Rcout << "Size " << L.size() << "\n";
	double tot=0;
	Rcpp::Rcout << "Do sum\n";
	for (int i=0;i<L.size();i++) {
		L[i]=exp(L[i]);
		tot=tot+L[i];
	}
	Rcpp::Rcout << tot << "\n";
	Rcpp::Rcout << "Do sort\n";
	sort(L.begin(),L.end());
	reverse(L.begin(),L.end());
	int found95=0;
	int found99=0;
	int found59=0;
	Rcpp::Rcout << "Do search  " << L.size() << "\n";
	for (int i=1;i<L.size();i++) {
		double unlog=log(L[i]);
		L[i]=L[i]+L[i-1];
		if (L[i]>tot*0.95&&found95==0) {
			found95=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << L[i] << "\n";
		}
		if (L[i]>0.99*tot&&found99==0) {
			found99=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << L[i] << "\n";
		}
		if (L[i]>0.99999*tot&&found59==0) {
			found59=1;
			Rcpp::Rcout << "Threshold " << unlog << " " << L[i] << "\n";
			break;
		}
	}
}

