#include "basicmodel.h"
#include "likelihoods.h"
#include "io.h"
#include <string>
#include "Rcpp.h"

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	//Infectious potential relative to time of symptoms
	//Derived from Ashcroft with Bonhoeffer
	p.pa=97.18750;
	p.pb=0.268908;
	p.po=25.625;

	//Time of symptoms relative to date of infection
	//Lognormal distribution cited by He et al.
	p.smu=1.434065;
	p.ssigma=0.6612;

	//The following are learnt from CUH data.
	//They characterise the distribution of time from positive test to date of reporting symptoms
	p.ucta=2.5932152095707406;
	p.uctb=3.7760060663975437;
	p.ucto=3.112080041460921;
	p.uct_mean=6.67992;
	//Can create an option to optimise over the unknown times of becoming symptomatic
	p.rate=0.0008;
	p.seq_noise=0.41369;
	p.threshold=0;
	p.thresholdns=0;
	p.max_n=10;
	p.min_qual=0.8;
	p.error=0;
    p.ali_file="";
    p.pat_file="";
    p.mov_file="";
	p.ward_file="";
	p.pat_delim='/';
	p.mov_delim='.';
	p.diagnostic=0;
	p.calc_thresholds=0;
	p.hcw_location_default=0.5714286;
	p.pat_location_default=1;
	p.ward_format_old=0;
	p.use_all_seqs=0;
	p.symptom_uncertainty_calc=0;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--evo_rate")==0) {
			x++;
			p.rate=atof(argv[x]);
		} else if (p_switch.compare("--seq_noise")==0) {
			x++;
			p.seq_noise=atof(argv[x]);
		} else if (p_switch.compare("--diag")==0) {
			x++;
			p.diagnostic=atoi(argv[x]);
		} else if (p_switch.compare("--ali_file")==0) {
			x++;
			p.ali_file=argv[x];
		} else if (p_switch.compare("--pat_file")==0) {
			x++;
			p.pat_file=argv[x];
		} else if (p_switch.compare("--mov_file")==0) {
			x++;
			p.mov_file=argv[x];
		} else if (p_switch.compare("--ward_file")==0) {/*New*/
			x++;
			p.ward_file=argv[x];
		} else if (p_switch.compare("--sub_file")==0) {
			x++;
			p.sub_file=argv[x];
		} else if (p_switch.compare("--threshold")==0) {
			x++;
			p.threshold=atoi(argv[x]);
		} else if (p_switch.compare("--hcw_loc_default")==0) {
			x++;
			p.hcw_location_default=atof(argv[x]);
		} else if (p_switch.compare("--pat_loc_default")==0) {
			x++;
			p.pat_location_default=atoi(argv[x]);
		} else if (p_switch.compare("--maxn")==0) {
			x++;
			p.max_n=atoi(argv[x]);
		} else if (p_switch.compare("--min_qual")==0) {
			x++;
			p.min_qual=atof(argv[x]);
		} else if (p_switch.compare("--calc_thresholds")==0) {
			x++;
			p.calc_thresholds=atoi(argv[x]);
		} else if (p_switch.compare("--use_all_seqs")==0) {
			x++;
			p.use_all_seqs=atoi(argv[x]);
		} else if (p_switch.compare("--ward_format_old")==0) {
			x++;
			p.ward_format_old=atoi(argv[x]);
        } else {
			Rcpp::Rcout << "Incorrect usage\n ";
			//			exit(1);
		}
		p_switch.clear();
		x++;
	}
	p.rate=(p.rate*29900)/365.25;
	SetThreshold(p);
}

void ReadPatFromCSVNoSeq(run_params& p, vector<pat>& pdata, const vector<int>& sdist_interval, const vector<double>& sdist_prob) {
	if (p.diagnostic==1) {
		Rcpp::Rcout << "Read data for individuals\n";
	}
	ifstream csv_file;
	csv_file.open(p.pat_file.c_str());
	string str;
	int i=-1;
	int spwarn=0;
	while (getline(csv_file,str)) {
		i++;
		if (i>0&&p.error==0) {
			pat pt;
			//Edit string to remove "
			RemovePunc(str);

			//Split by commas
			vector<string> subs;
			SplitCommas(str,subs);
			RemoveSpaces(p.pat_file.c_str(),i,spwarn,subs);

			pt.code=subs[0];
			pt.type=atoi(subs[3].c_str());
			if (pt.type==3) {
				pt.hcw=1;
			} else {
				pt.hcw=0;
			}

			//Sort out the dates
			int temp[] = {1};
			vector<int> d(temp,temp+sizeof(temp) / sizeof(int));
			for (int j=0;j<d.size();j++) {
				vector<int> dmy;
				MakeDMY(p,p.pat_file.c_str(),d[j],subs,p.pat_delim,dmy);
				int day=DatetoDay(dmy);
				//Rcpp::Rcout << "Day " << day << "\n";
				if (d[j]==1) {
					pt.time_s.most_likely=day;
				}
			}
			if (subs[2]=="1") {
				pt.time_s.time.push_back(pt.time_s.most_likely);
				pt.time_s.prob.push_back(1);
			} else {
				if (p.symptom_uncertainty_calc==0) {
					//Correct the symptom onset time by the mean of the difference between symptom and positive test times
					pt.time_s.time.push_back(pt.time_s.most_likely-floor(p.uct_mean+0.5));
					pt.time_s.most_likely=pt.time_s.most_likely-floor(p.uct_mean+0.5);
					pt.time_s.prob.push_back(1);
				} else {
					double max=-1;
					int max_index=0;
					for (int j=0;j<sdist_interval.size();j++) {
						pt.time_s.time.push_back(pt.time_s.most_likely+sdist_interval[j]);
						pt.time_s.prob.push_back(sdist_prob[j]);
						if (sdist_prob[j]>max) {
							max=sdist_prob[j];
							max_index=j;
						}
					}
					pt.time_s.most_likely=pt.time_s.most_likely-sdist_interval[max_index];
				}
			}
			pdata.push_back(pt);
		}
	}
}

void ReadPatFromCSV(run_params& p, vector<pat>& pdata, const vector<int>& sdist_interval, const vector<double>& sdist_prob) {
	if (p.diagnostic==1) {
		Rcpp::Rcout << "Read data for individuals\n";
	}
	ifstream csv_file;
	csv_file.open(p.pat_file.c_str());
	string str;
	int i=-1;
	int spwarn=0;
	while (getline(csv_file,str)) {
		if (p.diagnostic==1) {
			Rcpp::Rcout << "Line " << i << " " << str << "\n";
		}
		i++;
		//if (!(csv_file >> str)) break;
		if (i>0&&p.error==0) {
			pat pt;
			//Edit string to remove "
			RemovePunc(str);

			//Split by commas
			vector<string> subs;
			SplitCommas(str,subs);
			if (subs.size()!=7) {
				Rcpp::Rcout << "Error in input file: " << p.pat_file << " Wrong number of columns\n";
				p.error=1;
			}
			RemoveSpaces(p.pat_file.c_str(),i,spwarn,subs);
			
			pt.code=subs[0];
			pt.code_match=subs[4];
			pt.type=atoi(subs[3].c_str());
			if (pt.type==3) {
				pt.hcw=1;
			} else {
				pt.hcw=0;
			}

			//Sort out the dates
			int temp[] = {1,5};
			vector<int> d(temp,temp+sizeof(temp) / sizeof(int));
			for (int j=0;j<d.size();j++) {
				vector<int> dmy;
				MakeDMY(p,p.pat_file.c_str(),d[j],subs,p.pat_delim,dmy);
				int day=DatetoDay(dmy);
				if (d[j]==1) {
					pt.time_s.most_likely=day;
				}
				if (d[j]==5) {
					pt.time_seq=day;
				}
			}
			if (subs[2]=="1") {
				pt.time_s.time.push_back(pt.time_s.most_likely);
				pt.time_s.prob.push_back(1);
			} else {
				if (p.symptom_uncertainty_calc==0) {
					//Correct the symptom onset time by the mean of the difference between symptom and positive test times
					pt.time_s.time.push_back(pt.time_s.most_likely-floor(p.uct_mean+0.5));
					pt.time_s.most_likely=pt.time_s.most_likely-floor(p.uct_mean+0.5);
					pt.time_s.prob.push_back(1);
				} else {
					double max=-1;
					int max_index=0;
					for (int j=0;j<sdist_interval.size();j++) {
						pt.time_s.time.push_back(pt.time_s.most_likely+sdist_interval[j]);
						pt.time_s.prob.push_back(sdist_prob[j]);
						if (sdist_prob[j]>max) {
							max=sdist_prob[j];
							max_index=j;
						}
					}
					pt.time_s.most_likely=pt.time_s.most_likely-sdist_interval[max_index];
				}
			}

			//Rcpp::Rcout << str << "\n";
			//Rcpp::Rcout << pt.code << " " << pt.code_match << " " << pt.hcw << " " << pt.type << " " << pt.time_s.most_likely << " " << pt.time_seq << " " << pt.time_s_cert << "\n";
			pdata.push_back(pt);
		}
	}
}

void RemoveSpaces(string file, int i, int& pr, vector<string>& subs) {
	for (int j=0;j<subs.size();j++) {
		//Remove spaces
		vector<int> to_rem;
		for (int k=0;k<subs[j].size();k++){
			if (subs[j].compare(k,1," ")==0) {
				to_rem.push_back(k);
			}
		}
		sort(to_rem.begin(),to_rem.end());
		reverse(to_rem.begin(),to_rem.end());
		for (int k=0;k<to_rem.size();k++) {
			if (pr==0) {
				Rcpp::Rcout << "Warning: Removing spaces from line " << i << " of " << file << "\n";
				pr=1;
			}
			subs[j].erase(subs[j].begin()+to_rem[k]);
		}
	}
}

void ReadHCWMovFromCSV(run_params& p, vector<pat>& pdata) {
	ifstream csv_file;
	csv_file.open(p.mov_file.c_str());
	string str;
	vector<int> dates;
	int i=-1;
	int spwarn=0;
	while (getline(csv_file,str)) {
		i++;
		//Edit string to remove "
		RemovePunc(str);
		//Split at commas
		vector<string> subs;
		SplitCommas(str,subs);
		RemoveSpaces(p.mov_file.c_str(),i,spwarn,subs);
		if (i==0) {
			//Date information
			for (int j=2;j<subs.size();j++) {
				vector<int> dmy;
				MakeDMY(p,p.mov_file.c_str(),j,subs,p.mov_delim,dmy);
				int day=DatetoDay(dmy);
				//Rcpp::Rcout << subs[j] << " " << day << "\n";
				dates.push_back(day);
			}
		} else {
			//Put information into pdat location information
			string pat=subs[0];
			string w=subs[1];
			/*if (subs[1].compare("B")==0) {
				w="WARD_47";//N.B. Wards 47 to 50 all linked - code as 47
			}
			if (subs[1].compare("C")==0) {
				w="WARD_10";
			}
			if (subs[1].compare("A")==0) {
				w="WARD_19";
			}
			if (subs[1].compare("D")==0) {
				w="WARD_36";
			}
			if (subs[1].compare("E")==0) {
				w="WARD_52";
			}*/

			int index=-1;
			for (int j=0;j<pdata.size();j++) {
				if (pdata[j].code==pat) {
					index=j;
					//Rcpp::Rcout << pat << "\n";
					break;
				}
			}
			if (index!=-1) {
				for (int j=2;j<subs.size();j++) {
					if (subs[j].compare(0,1,"Y")==0) {
						loc l;
						l.ward=w;
						l.date=dates[j-2];
						l.prob=1;
						pdata[index].locat.push_back(l);
					}
				}
			}
		}
	}
}

void ReadWardMovFromCSV(run_params& p, vector<pat>& pdata) {
	cout << "Read ward file\n";
	ifstream csv_file;
	csv_file.open(p.ward_file.c_str());
	string str;
	vector<int> dates;
	int i=-1;
	int spwarn=0;
	while (getline(csv_file,str)) {
		i++;
		if (i>0) {
			//Rcpp::Rcout << "String here " << str << "\n";
			//Edit string to remove "
			RemovePunc(str);
			//Split at commas
			vector<string> subs;
			SplitCommas(str,subs);
			RemoveSpaces(p.ward_file.c_str(),i,spwarn,subs);
			//Look for matching patient
			string pat=subs[0];
			int index=-1;
			for (int j=0;j<pdata.size();j++) {
				if (pdata[j].code==pat) {
					index=j;
					//Rcpp::Rcout << pat << "\n";
					break;
				}
			}
			if (index!=-1) {
				for (int j=4;j<subs.size();j=j+3) {
					loc l;
					l.ward=subs[j];
					/*if (l.ward.compare("WARD_48")==0) {
						l.ward="WARD_47";
					}
					if (l.ward.compare("WARD_49")==0) {
						l.ward="WARD_47";
					}
					if (l.ward.compare("WARD_50")==0) {
						l.ward="WARD_47";
					}*/
					if (l.ward.size()>0) {
						vector<int> dmy;
						MakeDMY(p,p.ward_file.c_str(),j+1,subs,p.pat_delim,dmy);
						//Rcpp::Rcout << dmy[0] << " " << dmy[1] << " " << dmy[2] << "\n";
						int day1=DatetoDay(dmy);
						dmy.clear();
						if (p.error==0) {
							MakeDMY(p,p.ward_file.c_str(),j+2,subs,p.pat_delim,dmy);
							int day2=DatetoDay(dmy);
							for (int k=day1;k<=day2;k++) {
								l.date=k;
								l.prob=1;
								pdata[index].locat.push_back(l);
							}
						}
					}
				}
			}
		}
	}
}

void ReadWardMovFromCSVTemp(run_params& p, vector<pat>& pdata) {
	ifstream csv_file;
	csv_file.open(p.ward_file.c_str());
	string str;
	vector<int> dates;
	int i=-1;
	int spwarn=0;
	while (getline(csv_file,str)) {
		i++;
		if (i>0&&p.error==0) {
			//Rcpp::Rcout << str << "\n";
			/*if (p.diagnostic==1) {
				Rcpp::Rcout << "Ward_file string " << str << "\n";
			}*/
			//Edit string to remove "
			RemovePunc(str);
			//Split at commas
			vector<string> subs;
			SplitCommas(str,subs);
			RemoveSpaces(p.ward_file.c_str(),i,spwarn,subs);

			//Look for matching patient
			string pat=subs[0];
			int index=-1;
			for (int j=0;j<pdata.size();j++) {
				if (pdata[j].code==pat) {
					index=j;
					//Rcpp::Rcout << pat << "\n";
					break;
				}
			}
			if (p.diagnostic==1) {
				Rcpp::Rcout << subs.size() << " " << index << "\n";
			}

			if (index!=-1) {
				for (int j=1;j<subs.size();j=j+4) {
					loc l;
					l.ward=subs[j];
					/*if (l.ward.compare("WARD_48")==0) {
						l.ward="WARD_47";
					}
					if (l.ward.compare("WARD_49")==0) {
						l.ward="WARD_47";
					}
					if (l.ward.compare("WARD_50")==0) {
						l.ward="WARD_47";
					}*/
					if (l.ward.size()>0) {
						vector<int> dmy;
						MakeDMY(p,p.ward_file.c_str(),j+1,subs,p.pat_delim,dmy);
						int day1=DatetoDay(dmy);
						dmy.clear();
						if (p.error==0) {
							MakeDMY(p,p.ward_file.c_str(),j+3,subs,p.pat_delim,dmy);
							int day2=DatetoDay(dmy);
							if (p.diagnostic==1) {
								Rcpp::Rcout << day1 << " " << day2 << "\n";
							}
							for (int k=day1;k<=day2;k++) {
								l.date=k;
								l.prob=1;
								pdata[index].locat.push_back(l);
							}
						}
					}
				}
			}
		}
	}
}


void RemovePunc(string& str) {
	//Edit string to remove "
	vector<int> rem;
	for (int j=0;j<str.size();j++) {
		if (str[j]=='"') {
			rem.push_back(j);
		}
	}
	reverse(rem.begin(),rem.end());
	for (int j=0;j<rem.size();j++) {
		str.erase(str.begin()+rem[j]);
	}
}

void SplitCommas(const string str, vector<string>& subs) {
	stringstream ss(str);
	while (ss.good() ) {
		string sr;
		getline(ss,sr,',');
		subs.push_back(sr);
	}
}

void MakeDMY (run_params& p, string file, const int j, const vector<string>& subs, char delim, vector<int>& dmy) {
	stringstream sss(subs[j]);
	while (sss.good()&&p.error==0) {
		string sr;
		getline(sss,sr,delim);
		vector<string> months;
		months.push_back("Jan");
		months.push_back("Feb");
		months.push_back("Mar");
		months.push_back("Apr");
		months.push_back("May");
		months.push_back("Jun");
		months.push_back("Jul");
		months.push_back("Aug");
		months.push_back("Sep");
		months.push_back("Oct");
		months.push_back("Nov");
		months.push_back("Dec");
		int found=0;
		int add=-1;
		for (int i=0;i<months.size();i++) {
			if (sr.compare(months[i])==0) {
				add=i+1;
				found=1;
				break;
			}
		}
		if (found==1) {
			dmy.push_back(add);
		} else {
			dmy.push_back(atoi(sr.c_str()));
		}
	}
	if (p.error==0) {
		if (dmy.size()!=3) {
			DateError(file,subs[j],delim);
			p.error=1;
		}
		if (dmy[2]>50&&dmy[2]<100) {
			YearOOR(file,subs[j]);
			p.error=1;
		}
		if (dmy[2]>99&&dmy[2]<2000) {
			YearOOR(file,subs[j]);
			p.error=1;
		}
		if (dmy[2]>2050) {
			YearOOR(file,subs[j]);
			p.error=1;
		}
	}
}

void DateError (string file, string subsj, char delim) {
	Rcpp::Rcout << "Error reading in date\n" << subsj << " in " << file << "\n";
	Rcpp::Rcout << "Format must be DD" << delim << "MM" << delim << "YY or DD" << delim << "MM" << delim << "YYYY\n";
}

void YearOOR (string file, string subsj) {
	Rcpp::Rcout << "Error: Year out of range in entry\n" << subsj << " in " << file << "\n";
	Rcpp::Rcout << "Possible fix: check date is in DD/MM/YY or DD/MM/YYYY format \n";
}

int DatetoDay (vector<int>& dmy) {
	//Rcpp::Rcout << dmy[0] << " " << dmy[1] << " " << dmy[2] << "\n";
	int day=dmy[0];
	int lp=1;
	if (dmy[2]>1999) {
		dmy[2]=dmy[2]-2000;
	}
	//Year calculation
	int diff=dmy[2]-20;
	day=day+(diff*365);
	if (diff>0) {
		int leap=diff/4;
		day=day+leap+1;
	}
	if (diff%4!=3) {
		lp=0;
	}
	//Month calculation
	if (dmy[1]>1) {
		day=day+31;
	}
	if (dmy[1]>2) {
		day=day+28+lp;
	}
	if (dmy[1]>3) {
		day=day+31;
	}
	if (dmy[1]>4) {
		day=day+30;
	}
	if (dmy[1]>5) {
		day=day+31;
	}
	if (dmy[1]>6) {
		day=day+30;
	}
	if (dmy[1]>7) {
		day=day+31;
	}
	if (dmy[1]>8) {
		day=day+31;
	}
	if (dmy[1]>9) {
		day=day+30;
	}
	if (dmy[1]>10) {
		day=day+31;
	}
	if (dmy[1]>11) {
		day=day+30;
	}
	if (dmy[1]>12) {
		day=day+31;
	}
	//Rcpp::Rcout << "Day " << day << "\n";
	return day;
}

void EditHCWMovData (vector<pat>& pdat) {
	//Adds a 12-hour window of uncertainty in the location of HCWs (not patients)
	//Represents e.g. transmission via touching equipment, night shift timings, etc.
	//Assign a probability of 0.5 to the days before and after a known presence
	for (int i=0;i<pdat.size();i++) {
		vector<loc> new_loc;
		for (int j=0;j<pdat[i].locat.size();j++) {
			loc l=pdat[i].locat[j];
			int minus=1;
			int plus=1;
			if (j>0&&pdat[i].locat[j-1].date==l.date-1&&pdat[i].locat[j-1].ward==l.ward&&pdat[i].locat[j-1].prob==1) {
				minus=0;
			}
			if (j<pdat[i].locat.size()-1&&pdat[i].locat[j+1].date==l.date+1&&pdat[i].locat[j+1].ward==l.ward&&pdat[i].locat[j+1].prob==1) {
				plus=0;
			}
			if (minus==1) {
				loc lm=l;
				lm.date--;
				lm.prob=0.5;
				new_loc.push_back(lm);
			}
			if (plus==1) {
				loc lp=l;
				lp.date++;
				lp.prob=0.5;
				new_loc.push_back(lp);
			}
		}
		for (int j=0;j<new_loc.size();j++) {
			pdat[i].locat.push_back(new_loc[j]);
		}
	}
}

void ReadFastaAli (run_params p, vector<string>& names, vector<string>& seqs) {
    ifstream ali_file;
    ali_file.open(p.ali_file.c_str());
	string seq;
	string str;
	for (int i=0;i<1000000;i++) {
		if (!(ali_file >> str)) break;
		if (str.at(0)=='>') {
			names.push_back(str);
			if (p.diagnostic==1) {
				Rcpp::Rcout << "Read seqname " << str << "\n";
			}
			if (seq.size()>0) {
				seqs.push_back(seq);
				seq.clear();
			}
		} else {
			seq=seq+str;
		}
	}
	if (seq.size()>0) {
		seqs.push_back(seq);
	}
}

void CheckBaseCase (vector<string>& seqs) {
	Rcpp::Rcout << "Check case of sequence data\n";
	for (int i=0;i<seqs.size();i++) {
		//Rcpp::Rcout << "Check sequence " << i << " " << seqs[i].size() << "\n";
		for (int j=0;j<seqs[i].size();j++) {
			//Rcpp::Rcout << "Check character " << i << " " << j << " " << seqs[i][j] << "\n";
			if (seqs[i].compare(j,1,"a")==0) {
				seqs[i][j]='A';
			}
			if (seqs[i].compare(j,1,"c")==0) {
				seqs[i][j]='C';
			}
			if (seqs[i].compare(j,1,"g")==0) {
				seqs[i][j]='G';
			}
			if (seqs[i].compare(j,1,"t")==0) {
				seqs[i][j]='T';
			}
			if (seqs[i].compare(j,1,"n")==0) {
				seqs[i][j]='N';
			}
		}
	}
	Rcpp::Rcout << "CheckBaseCase complete\n";
}

void ReadSubsets (run_params p, vector< vector<int> >& subsets) {
	ifstream str_file;
	//Rcpp::Rcout << p.sub_file << "\n";
    str_file.open(p.sub_file.c_str());
	string str;
	for (int i=0;i<100000;i++) {
		if (!(str_file >> str)) break;
		//Rcpp::Rcout << str << "\n";
		vector<string> subs;
		SplitCommas(str,subs);
		vector<int> vals;
		for (int j=0;j<subs.size();j++) {
			vals.push_back(atoi(subs[j].c_str()));
			//Rcpp::Rcout << atoi(subs[j].c_str()) << "\n";
		}
		subsets.push_back(vals);
	}
}

