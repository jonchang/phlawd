#ifndef FASTA_UTIL_H_
#define FASTA_UTIL_H_

#define FULL_METADATA -987
#define TAXON_NAMES -876
#define NCBI_TAXON_IDS -765

#include <string>

using namespace std;

class FastaUtil{
private:
	bool readFile(string &, vector<Sequence> &, bool is_aligned, bool is_user_fasta);

public:
	FastaUtil();
	bool read_aligned_fasta_into(vector<Sequence> & seqs, string & filen);
	bool read_unaligned_fasta_into(vector<Sequence> & seqs, string & filen);
	bool read_user_fasta_into(vector<Sequence> & seqs, string & filen);

	bool writeFileFromVector(string filen, vector<Sequence> & seqs, int label_format = FULL_METADATA);
};

#endif
