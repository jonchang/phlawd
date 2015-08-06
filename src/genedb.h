#ifndef _GENEDB_H
#define _GENEDB_H

#include <string>
#include <vector>
#include <map>

#include "utils.h"
#include "sequence.h"
#include "fasta_util.h"

using namespace std;

class GeneDB{
private:
    string name;//database filename

	void begin_transaction(sqlite3 * db);
	void commit_transaction(sqlite3 * db);
	void execute_query(sqlite3 * db, string & query); // this one is just a wrapper for the const char implementation
	void execute_query(sqlite3 * db, const char * query);
	void execute_simple_transaction(sqlite3 * db, string & query);
	int insert_sequence(sqlite3 * db, Sequence & seq);
	int insert_alignment(sqlite3 * db, string & alignment_name, int updated);
	int insert_profile_alignment(sqlite3 * db, string & profile_name, int child1_id, int child2_id);
	int insert_sequence_alignment_map(sqlite3 * db, int sequence_id, int alignment_id, string & sequence);
	int insert_sequence_profile_map(sqlite3 * db, int sequence_id, int profile_id, string & sequence);
	void load_seqs_from_alignment_into(vector<Sequence> & seqs, int alignment_id, bool is_profile, bool want_aligned);


public:
    GeneDB();
    GeneDB(string name);

    // adding and removing db records
    // these are low level queries primarily used within the GeneDB class, should perhaps be made private
    void add_sequences(vector<Sequence> * seqs);
//    int add_original_alignment(string & filen, vector<Sequence> * dbseqs, vector<Sequence> * userseqs);
    int add_original_alignment(string & filen, vector<Sequence> * seqs);
    int add_profile_alignment(string & profile_name, int child1_id, int child2_id, vector<Sequence> & sequences);
    int add_profile_alignment(string & profile_name, vector<Sequence> & sequences);
    int add_empty_intermediate_profile(int child1_id, int child2_id); // TODO: just for SQLiteProfiler, should be deprecated
    void remove_original_alignment_by_id(int alignid);
    void remove_original_alignment_by_name(string alignname);
    void remove_profile_alignments();

    // retrieving information and data from the db
    string get_db_name();
    string get_original_alignment_name_for_db_id(int alignment_id);
    void get_align_seq_unaligned_fully_initialized(string alignname, vector<Sequence> & seqs);
    void get_first_profile_alignments(vector<string> & names);
    int get_original_alignment_id_by_name(string & alignname);
    Sequence create_full_unaligned_sequence_obj_for_original_sequence_db_id(int sqlite_id);
    void load_all_original_sequences_into(vector<Sequence> & seqs);
    void load_orig_alignment_labels_into(vector<string> & names);
    void load_orig_alignment_db_ids_into(vector<int> & nums);
    void load_first_profile_ids_into(vector<int>&);

    void load_aligned_seqs_from_profile_alignment_into(vector<Sequence> & seqs, int profile_id);
    void load_unaligned_seqs_from_profile_alignment_into(vector<Sequence> & seqs, int profile_id);
    void load_aligned_seqs_from_original_alignment_into(vector<Sequence> & seqs, int alignment_id);
    void load_unaligned_seqs_from_original_alignment_into(vector<Sequence> & seqs, int alignment_id);

    void load_original_alignment_info_into(map<int, string> & the_map);
    void load_alignments_to_update_info_into(map<int, string> & the_map);

    // miscellaneous (mostly higher-level) db modification tasks
    void initialize(bool overwrite);
    void toggle_alignment_update(int alignid);
    void associate_sequence_with_alignment(int alignid,Sequence & inseq);
    void update_align_seqs(int alignid,vector<Sequence> & seqs);
    void update_profile_align_seqs(int alignid, vector<Sequence> & seqs);
    void migrate_original_alignments_and_load_info_into(map<int,string> & profile_id_name_map);
    void migrate_alignments_to_update_and_load_info_into(map<int,string> & profile_id_name_map, vector<int> & updatedprofsnums);

    // not yet addressed
/*    void copy_alignments_to_first_profiles_updated(map<int, string> & profile_id_name_map,vector<int> & updatedprofsnums); */
    int get_deepest_profile_for_alignment(int alignid);
    void add_sequences_for_profile_alignment(int profilealignid,vector<Sequence> & seqs);
    void write_profile_alignment_to_file(int alignid, string filename, int format = FULL_METADATA);
//    void write_profile_alignment_with_names_to_file(int alignid, string filename,bool ncbi);
    void get_updated_profs_names_delete_old(vector<string> & updatedprofs,vector<int> & updatedprofsnums,vector<int> & notupdatedprofnums);
    void toggle_updated_all_off();
};

#endif
