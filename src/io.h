#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

string get_extdata(string fname);

void GetOptions (run_params& p, int argc, const char **argv);
void GetThresholds (vector< vector<double> >& thresholds95, vector< vector<double> >& thresholds99, double& t95NS, double& t99NS, int& error);
void ReadPatFromCSV(run_params& p, vector<pat>& pdata, const vector<int>& sdist_interval, const vector<double>& sdist_prob);
void ReadPatFromCSVNoSeq(run_params& p, vector<pat>& pdata, const vector<int>& sdist_interval, const vector<double>& sdist_prob);
void RemoveSpaces(string file, int i, int& pr, vector<string>& subs);
void ReadHCWMovFromCSV(run_params& p, vector<pat>& pdata);
void ReadWardMovFromCSV(run_params& p, vector<pat>& pdata);
void ReadWardMovFromCSVTemp(run_params& p, vector<pat>& pdata);
void EditHCWMovData (vector<pat>& pdat);

void RemovePunc(string& str);
void SplitCommas(const string str, vector<string>& subs);
void MakeDMY (run_params& p, string file, const int j, const vector<string>& subs, char delim, vector<int>& dmy);
void DateError (string file, string subsj, char delim);
void YearOOR (string file, string subsj);

int DatetoDay (vector<int>& dmy);
void ReadFastaAli (run_params p, vector<string>& names, vector<string>& seqs);
void CheckBaseCase (vector<string>& seqs);
void ReadSymptomData (vector<pat>& pdat);
void ReadSeqTimes(vector<pat>& pdat);
void ReadLocationData (vector<pat>& pdat);
void ReadCommunityDistances (vector<int>& comm_dist);
void ReadSubsets (run_params p, vector< vector<int> >& subsets);
