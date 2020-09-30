#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void FindVariants (vector<sparseseq>& variants, string& consensus, vector<pat>& pdat);
void FindVariantsNoSeq (vector<sparseseq>& variants, vector<pat>& pdat);
void ListAllVariantPositions (const vector<sparseseq>& variants, vector<allele>& allvar);
void FindAmbiguousVarPositions (const vector<allele>& allvar, vector<pat>& pdat, vector<int>& nloc_count);
bool compare_allele(allele v1, allele v2);




