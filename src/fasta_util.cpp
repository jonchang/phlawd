#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

#include "sequence.h"
#include "fasta_util.h"
#include "utils.h"

FastaUtil::FastaUtil() {
}

struct upper {
	int operator()(int c) {
		return std::toupper((unsigned char) c);
	}
};

bool FastaUtil::read_aligned_fasta_into(vector<Sequence> & seqs, string & filen) {
	bool is_aligned = true;
	bool is_user_fasta = false;
	return readFile(filen, seqs, is_aligned, is_user_fasta);
}

bool FastaUtil::read_unaligned_fasta_into(vector<Sequence> & seqs, string & filen) {
	bool is_aligned = false;
	bool is_user_fasta = false;
	return readFile(filen, seqs, is_aligned, is_user_fasta);
}

bool FastaUtil::read_user_fasta_into(vector<Sequence> & seqs, string & filen) {
	bool is_aligned = false;
	bool is_user_fasta = true;
	return readFile(filen, seqs, is_aligned, is_user_fasta);
}

// TODO: return false if not a fasta
bool FastaUtil::readFile(string & filen, vector<Sequence> & seqs, bool is_aligned, bool is_user_fasta) {
	string tline;
	ifstream infile(filen.c_str());
	bool first = true;
	Sequence cur;
	string curseq;

	while (getline(infile, tline)) {

		// skip blank lines
		TrimSpaces(tline);
		if (tline.size() < 1) {
			continue;
		}

		if (tline.substr(0, 1) == ">") { // this line defines a new sequence
			string label_ = tline.substr(1, tline.size() - 1);

			if (first == true) {
				// first seq in file, don't try to record the last one (it's not there)
				first = false;

			} else { // not the first seq in the file, so record the last one

				// make the last seq uppercase, add it to the last seq object
				std::transform(curseq.begin(), curseq.end(), curseq.begin(), upper());
				if (is_aligned)
					cur.set_aligned_sequence(curseq);
				else
					cur.set_unaligned_sequence(curseq);

				// put the last (now complete) seq object in the stack
				seqs.push_back(cur);
			}

			// start a new seq object
			cur = Sequence(label_, is_user_fasta);
			curseq = "";

		} else {
			// keep adding sequence lines till we get to the next sequence
			curseq += tline;
		}
	}

	// finish the last seq object before we stop
	std::transform(curseq.begin(), curseq.end(), curseq.begin(), upper());
	if (is_aligned)
		cur.set_aligned_sequence(curseq);
	else
		cur.set_unaligned_sequence(curseq);

	// add the last seq to the stack
	seqs.push_back(cur);

	// done
	infile.close();
	return true;
}

/* this is just bare bones, write a vector of sequences to a file.
 *
 * the format is for the labels. it can be:
 *
 * 0 - taxon ids only (user seqs have taxon names printed)
 * 1 - taxon names
 * 2 - full label metadata see the Sequence class for format -- DEFAULT
 */
bool FastaUtil::writeFileFromVector(string filename, vector<Sequence> & seqs, int label_format) {

	ofstream outfile;
	outfile.open(filename.c_str(), ios::out);
	for (unsigned int i = 0; i < seqs.size(); i++) {

		/*
		// write the sequence label
		outfile << ">";
		if (seqs[i].is_user_seq())
			outfile << "user_";
		outfile << seqs[i].get_label();
		*/

		// write the sequence label
		outfile << ">";
		if (label_format == FULL_METADATA)
			outfile << seqs[i].get_label();

		else if (label_format == TAXON_NAMES)
			outfile << seqs[i].get_taxon_name();

		else if (label_format == NCBI_TAXON_IDS)
			if (seqs[i].is_user_seq())
				outfile << "user_" + seqs[i].get_taxon_name();
			else
				outfile << seqs[i].get_ncbi_tax_id();

		else {
			cerr << "invalid output format specified. phlawd will exit";
			exit(1);
		}

		// write the sequence
		outfile << "\n";
		outfile << seqs[i].get_sequence();
		outfile << "\n";
	}
	outfile.close();
	return true;
}
