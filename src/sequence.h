#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#define SEQ_FIELD_SEP "|"

#include <string>

using namespace std;

class Sequence{
private:

/*	string label; /* the use of the label field currently differs between sequences of ncbi
	 * origin and user-inputted sequences.
	 *
	 * for ncbi sequences, we (should in general) store the ncbi taxon id in the label (this
	 * may change, it would perhaps be better to simply store the taxon name or something else
	 * entirely), while for user sequences, we store the identifier string from the user's
	 * fasta file, since there may not be a corresponding ncbi taxon for a given user sequence.
	 *
	 * one important use of the label field is for outputting fasta files; it forms the
	 * identifier string that follows the '>' character indicating the start of each sequence.
	 *
	 * note that when outputting user sequences to fasta files, we must differentiate their
	 * labels from potential ncbi taxon names in some way that is guaranteed not to generate
	 * duplicate taxon names. the current practice is to append 'user_' to the beginning of
	 * each user sequence label (see.
	 */

	string ncbi_gi;				// the ncbi accession number for this sequence
    string ncbi_tax_id;			// the ncbi taxon id of the taxon this sequence is assigned to
    string unaligned_seq;		// the unaligned sequence (i.e. no gaps!)
    string aligned_seq;			// the aligned sequence
    string tax_name;			// taxon name (in general we use the ncbi taxonomy to set this)
    int sqlite_id;				// the sqlite3 database id of this sequence in the phlawd db
    bool b_is_aligned;			// whether this sequence has been aligned, if this is true we return the aligned seq by default
    bool b_is_user_seq;			// whether this is a user sequence
    bool b_is_exemplar;			// whether this has been renamed to exemplify a higher taxon (see find_best_exemplar_for_higher_taxon in SQLiteController)

    bool b_is_original_seq;		// set to true for original sequences read in from fasta (e.g. original user sequences)
    string original_label;		// only used for original sequences, since they do not contain phlawd metadata in the label

    string comment;				// whatever else we want... not sure how this is used though

    string reverse(string);		// make the reverse complement of this sequence

    void set_is_original_seq(bool _is_original_seq);
    void set_original_label(string _original_label);
    void parse_label(string label, bool is_user_fasta);

public:
    Sequence();
    Sequence(string label, bool is_user_fasta);
/*    Sequence(int,string);
    Sequence(string,string,bool);
    Sequence(string,string); */
    bool is_aligned();
    bool is_user_seq();
    bool is_original_seq();
    bool is_exemplar();
    int get_sqlite_id();
    string get_sequence();
    string get_aligned_seq();
    string get_unaligned_seq();
    string get_ncbi_tax_id();
    string get_ncbi_gi_number();
    string get_label();
    string get_taxon_name();
    string get_comment();
    string reverse_complement();

    Sequence clone();

    void set_unaligned_sequence(string seq);
//    void set_label(string label);
    void set_taxon_name(string name);
    void set_comment(string comment);
    void set_is_aligned(bool al);
    void set_is_user_seq(bool user_seq);
    void set_is_exemplar(bool _is_exemplar);
    void perm_reverse_complement();
    void set_sqlite_id(int iid);
    void set_aligned_sequence(string inseq);
    void set_ncbi_tax_id(string _tid);
    void set_ncbi_gi_number(string _tid);
    bool operator==(const Sequence &other) const;
};
#endif
