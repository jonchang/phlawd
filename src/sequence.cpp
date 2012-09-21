/*
 * This is an extremely stripped down sequence class that is just meant to be
 * transparent and lightweight. As functionality increases, so will the 
 * complexity of the class.
 *
 */

#include <string>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include "utils.h"

using namespace std;

#include "sequence.h"

template<class T>
inline std::string to_string(const T & t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}
/* create an empty sequence */
Sequence::Sequence() :
		sqlite_id(-1),
		original_label(""),
		tax_name(""),
		ncbi_gi(""),
		ncbi_tax_id(""),
		aligned_seq(""),
		unaligned_seq(""),
		b_is_aligned(false),
		b_is_user_seq(false),
		b_is_original_seq(false),
		b_is_exemplar(false),
		comment("") {
}

/* accepts a phlawd sequence label string that is parsed
 * to populate the new sequence's metadata fields */
Sequence::Sequence(string label, bool is_user_fasta = false) :
		sqlite_id(-1),
		original_label(""),
		tax_name(""),
		ncbi_gi(""),
		ncbi_tax_id(""),
		aligned_seq(""),
		unaligned_seq(""),
		b_is_aligned(false),
		b_is_user_seq(false),
		b_is_original_seq(false),
		b_is_exemplar(false),
		comment("") {

	// populate the metadata fields with the info from the label
	parse_label(label, is_user_fasta);

}

/* returns a new sequence object containing identical
 * identical information to this sequence */
Sequence Sequence::clone() {

	Sequence cloned_seq = Sequence();
	cloned_seq.set_sqlite_id(sqlite_id);
	cloned_seq.set_original_label(original_label);
	cloned_seq.set_taxon_name(tax_name);
	cloned_seq.set_ncbi_gi_number(ncbi_gi);
	cloned_seq.set_ncbi_tax_id(ncbi_tax_id);
	cloned_seq.set_aligned_sequence(aligned_seq);
	cloned_seq.set_unaligned_sequence(unaligned_seq);
	cloned_seq.set_is_aligned(b_is_aligned);
	cloned_seq.set_is_user_seq(b_is_user_seq);
	cloned_seq.set_is_original_seq(b_is_original_seq);
	cloned_seq.set_is_exemplar(b_is_exemplar);
	cloned_seq.set_comment(comment);

	return(cloned_seq);
}

/*
 * sequence constructors have been deprecated to make sure the assignment of sequence fields is explicit
 *
 *
Sequence::Sequence(int _sqlite_id, string _al_seq) :
		sqlite_id(_sqlite_id), aligned_seq(_al_seq), aligned(), name(""), comment(""), ncbi_gi("0"), ncbi_tax("0") {
}
Sequence::Sequence(string _id, string _seq, bool _aligned) :
		id(_id), unaligned_seq(_seq), aligned(_aligned), name(""), comment(""), ncbi_gi("0"), ncbi_tax("0") {
}
Sequence::Sequence(string _id, string _seq) :
		id(_id), unaligned_seq(_seq), aligned(false), name(""), comment(""), ncbi_gi("0"), ncbi_tax("0") {
} */

////////////////////////////////////////////////////////////////////
//
//		private methods
//
////////////////////////////////////////////////////////////////////

void Sequence::set_is_original_seq(bool _is_original_seq) {
	b_is_original_seq = _is_original_seq;
}

void Sequence::set_original_label(string _original_label) {
	original_label = _original_label;
}

void Sequence::parse_label(string label, bool is_user_fasta) {

	/* populate the sequence metadata fields with information from the label
	 *
	 * format: (also used by get_label to create sequence labels):
	 *
	 * source | tax_name | ncbi_tax_id | gi | sqlite_id | is_aligned
	 */

	vector<string> toks = split(label, SEQ_FIELD_SEP);

	// case 1: when reading in the user sequences, we treat the entire label as the taxon name
	if (is_user_fasta) {
		set_taxon_name(label);
		set_is_user_seq(true);

	// case 2: when reading in an alignment with phlawd-formatted labels
	} else if (toks[0].compare("user") == 0 || toks[0].compare("ncbi") == 0) {

		if (toks[0].compare("user") == 0)
			set_is_user_seq(true);
		else
			set_is_user_seq(false);

		// strings are read straight in
		set_taxon_name(toks[1]);
		set_ncbi_tax_id(toks[2]);
		set_ncbi_gi_number(toks[3]);

		// have to convert numeric values
		set_sqlite_id(atoi(toks[4].c_str()));
		set_is_aligned(!!atoi(toks[5].c_str())); // expecting boolean ('1' or '0') for !!

	// case 3: a label that doesn't appear to be generated by phlawd (but we are not reading as a user seq)
	// ... this is not actually implemented anywhere yet.
	} else {

		// don't attempt to parse the label
		set_is_original_seq(true);
		set_original_label(label);
	}
}

string Sequence::reverse(string charin) {

	/* this should eventually be generalized for both nucleotide and
	 * protein, but for now it is just nucleotide for simplicity.
	 */

	string ret;
	if (charin == "-")
		ret = "-";
	if (charin == " ")
		ret = " ";
	if (charin == "A" || charin == "a") {
		ret = "T";
	} else if (charin == "T" || charin == "t") {
		ret = "A";
	} else if (charin == "C" || charin == "c") {
		ret = "G";
	} else if (charin == "G" || charin == "g") {
		ret = "C";
	} else if (charin == "U" || charin == "u") {
		ret = "A";
	} else if (charin == "m" || charin == "M") {
		ret = "K";
	} else if (charin == "r" || charin == "R") {
		ret = "Y";
	} else if (charin == "y" || charin == "Y") {
		ret = "R";
	} else if (charin == "k" || charin == "K") {
		ret = "M";
	} else if (charin == "v" || charin == "V") {
		ret = "B";
	} else if (charin == "h" || charin == "H") {
		ret = "D";
	} else if (charin == "d" || charin == "D") {
		ret = "H";
	} else if (charin == "b" || charin == "B") {
		ret = "V";
	} else if (charin == "n" || charin == "N" || charin == "x" || charin == "X") {
		ret = "N";
	}
	return ret;
}

////////////////////////////////////////////////////////////////////
//
//		getter methods
//
////////////////////////////////////////////////////////////////////

bool Sequence::is_aligned() {
	return b_is_aligned;
}

bool Sequence::is_user_seq() {
	return b_is_user_seq;
}

bool Sequence::is_original_seq() {
	return b_is_original_seq;
}

bool Sequence::is_exemplar() {
	return b_is_exemplar;
}

int Sequence::get_aligned_length() {
	return get_aligned_seq().size();
}

int Sequence::get_unaligned_length() {
	return get_unaligned_seq().size();
}

string Sequence::get_sequence() {

	if (is_aligned())
		return get_aligned_seq();
	else
		return get_unaligned_seq();

}

string Sequence::get_source() {
	if (is_user_seq())
		return "user";
	else
		return "ncbi";
}

string Sequence::get_aligned_seq() {
	if (is_aligned() == false) {
		cerr << "invalid request for an aligned sequence from an unaligned sequence object" << endl;
		exit(1);

	} else
		return aligned_seq;
}

string Sequence::get_unaligned_seq() {
	string seq = unaligned_seq;
	if (seq.compare((string)"") == 0) {
		cerr << "invalid request for an unaligned sequence from a sequence object that does not contain one" << endl;
		exit(1);

	} else
		return unaligned_seq;
}

/*
string Sequence::get_label() {
	return label;
} */

string Sequence::get_label() {

	/* generates a label string for this sequence and returns it.
	 *
	 * format (also used by parse_label to populate sequence fields):
	 *
	 * source | tax_name | ncbi_tax_id | gi | sqlite_id | is_aligned
	 */

	string label = "";

	if (is_original_seq()) {
		label = original_label;

	} else { // this is a phlawd-generated sequence

		if (is_user_seq())
			label += "user";
		else
			label += "ncbi";

		label += SEQ_FIELD_SEP + get_taxon_name();
		label += SEQ_FIELD_SEP + get_ncbi_tax_id();
		label += SEQ_FIELD_SEP + get_ncbi_gi_number();
		label += SEQ_FIELD_SEP + to_string(get_sqlite_id());
		label += SEQ_FIELD_SEP + to_string(is_aligned());
	}

	return label;
}

string Sequence::get_taxon_name() {
	return tax_name;
}

string Sequence::get_comment() {
	return comment;
}

string Sequence::get_ncbi_tax_id() {
	return ncbi_tax_id;
}

string Sequence::get_ncbi_gi_number() {
	return ncbi_gi;
}

int Sequence::get_sqlite_id() {
	return sqlite_id;
}

////////////////////////////////////////////////////////////////////
//
//		setter methods
//
////////////////////////////////////////////////////////////////////

void Sequence::set_unaligned_sequence(string _seq) {
	unaligned_seq = _seq;
}

void Sequence::set_aligned_sequence(string inseq) {
	aligned_seq = inseq;
	set_is_aligned(true);
}

/*
void Sequence::set_label(string _label) {
	label = _label;
} */

void Sequence::set_is_aligned(bool _is_aligned) {
	b_is_aligned = _is_aligned;
}

void Sequence::set_is_user_seq(bool _is_user_seq) {
	b_is_user_seq = _is_user_seq;
}

void Sequence::set_is_exemplar(bool _is_exemplar) {
	b_is_exemplar = _is_exemplar;
}

void Sequence::set_taxon_name(string _tax_name) {
	tax_name = _tax_name;
}

void Sequence::set_comment(string _comment) {
	comment = _comment;
}

void Sequence::set_ncbi_tax_id(string _ncbi_tax_id) {
	ncbi_tax_id = _ncbi_tax_id;
}

void Sequence::set_ncbi_gi_number(string _ncbi_gi) {
	ncbi_gi = _ncbi_gi;
}

void Sequence::set_sqlite_id(int iid) {
	sqlite_id = iid;
}

////////////////////////////////////////////////////////////////////
//
//		other public methods
//
////////////////////////////////////////////////////////////////////


string Sequence::reverse_complement() {
	string rcomp = get_sequence();
	string seq = get_sequence();

	for (unsigned int i = 0; i < rcomp.size(); i++) {
		rcomp.replace(i, 1, reverse(seq.substr(seq.size() - i - 1, 1)));
	}
	return rcomp;
}

void Sequence::perm_reverse_complement() {
	string rcomp = get_sequence();
	string seq = get_sequence();

	for (unsigned int i = 0; i < rcomp.size(); i++) {
		rcomp.replace(i, 1, reverse(seq.substr(seq.size() - i - 1, 1)));
	}
	seq = rcomp;
}

bool Sequence::operator==(const Sequence &other) const {
	bool ret = true;
	if (other.ncbi_gi != ncbi_gi) {
		ret = false;
		return ret;
	}
	if (other.ncbi_tax_id != ncbi_tax_id) {
		ret = false;
		return ret;
	}
	if (other.tax_name != tax_name) {
		ret = false;
		return ret;
	}
	return ret;
}
