/*
 * PHLAWD: Phylogeny assembly with databases
 * Copyright (C) 2010  Stephen A. Smith
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

/*
 *  SQLiteConstructor.h
 */

#ifndef _SQLITECONSTRUCTOR_H_
#define _SQLITECONSTRUCTOR_H_

#include <string>
#include <vector>
#include <map>

using namespace std;

#include "libsqlitewrapped.h"

#include "sequence.h"

#include "tree.h"
#include "node.h"
#include "genedb.h"

struct seqRecordTuple;
typedef struct seqRecordTuple seqRecordTuple;

class SQLiteConstructor {
private:
    ofstream logfile;
    fstream gifile;
    fstream ufafile;
    string clade_name;
    vector <string> search;
    bool searchliteral;
    string gene_name;
    string gene_db_name;
    string genefoldername;
    GeneDB gene_db;
    double mad_cutoff;
    double coverage;
    double identity;
    string db;
    int port;
    string listfilename;
    string excludefilename;
    string exclude_gi_filename;
    string include_gi_filename;
    bool onlynamesfromfile;
    bool excludenamesfromfile;
    bool includegifromfile;
    bool excludegifromfile;
    bool containshigher;
    bool containswild;
    bool containshigherex;
    bool containswildex;
    bool useITS;
    int numthreads;
    bool automated;
    bool updateDB;
    bool updateFILE;
    string updatef;
    bool ncbi_saturation;
    bool usertree;
    string usertreefile;
    Tree * userguidetree;
    bool userfasta;
    string userfastafile;
    bool userskipsearchdb;
    bool skipdbcheck;
    bool justseqquery;
    bool assignleftovers;
    int shrinkablethreshold;
    int main_left; // left value for clade id
    int main_right; // right value for clade id
    map<Sequence*,Node*> user_fasta_node_map; 
    vector<Sequence> * user_seqs;
    vector<Sequence> * known_seqs;
    vector<Sequence> use_only_names_from_file(vector<Sequence>& seqs);

    Sequence find_best_exemplar_for_higher_taxon(string taxon_id,vector<Sequence>& seqs); // TODO: needs a better name
    void load_sequences_filtered_by_ncbi_taxon_id_into(vector<Sequence> & filtered_seqs, vector<Sequence> & seqs_to_filter, vector<string> & ncbi_taxon_ids);

    vector<Sequence> exclude_names_from_file(vector<Sequence>& seqs);
    vector<Sequence> exclude_gis_from_file(vector<Sequence>& seqs);
    vector<Sequence> include_gis_from_file(vector<Sequence>& seqs);

    void load_info_for_all_seqs_matching_search_into(vector<seqRecordTuple> & found_seqs);
    vector<Sequence> make_seqs_from_seq_tuples_for_taxon(int search_clade_id, vector<seqRecordTuple> & results);

    void get_best_hits_openmp_SWPS3(vector<Sequence> & seqs,  vector<Sequence> * keep_seqs);
    void get_best_hits_openmp_SWPS3_justquery(vector<Sequence> & seqs,  vector<Sequence> * keep_seqs);
    double get_usertree_keepseq_overlap(vector<Sequence> * keep_seqs);
    void remove_duplicates_SWPS3(vector<Sequence> * keep_seqs);
    void reduce_genomes(vector<Sequence> * keep_seqs);
    void clean_shrunken_genomes();

    void load_info_on_valid_immediate_children_of_taxon_id_into(string parent_taxon_id, vector<string> & child_ids, vector<string> & child_tax_names);
    void find_db_child_seqs_of_ncbi_taxon_id(string parent_taxon_id, vector<Sequence> * seqs_to_search, vector<Sequence> * found_seqs);
    vector<string> get_valid_ncbi_child_taxon_ids_for_parent_id(string name_id);

    void remove_seq_from_vector_by_ncbi_id(vector<Sequence> * inseqs, string ncbi_id);
    void remove_seq_from_vector_by_taxon_name(vector<Sequence> * inseqs, string tname);

    void get_seqs_for_names_user(string name_id, vector<Sequence> * seqs);
    void get_seqs_for_nodes(Node * node, vector<Sequence> * seqs, vector<Sequence> * temp_seqs);
    void get_seqs_for_user_nodes(Node * node, vector<Sequence> * seqs);
    vector<string> get_final_children_node(Node * node);
    vector<string> get_final_children_node_hier(Node * node);

    void make_mafft_multiple_alignment(vector<Sequence> * inseqs);
    void make_mafft_multiple_alignment(vector<Sequence> * inseqs, vector<Sequence> * inseqs2);


    double calculate_MAD_quicktree();
    double calculate_MAD_quicktree_sample(vector<Sequence> * inseqs, vector<Sequence> * inuserseqs);
    void saturation_tests(vector<string> name_ids, vector<string> names, vector<Sequence> * keep_seqs);
    int get_single_to_group_seq_score(Sequence & inseq,vector<Sequence> & ginseqs);
    void write_gi_numbers(vector<Sequence> *);
    void write_user_numbers();
    void update_seqs_using_last_alignment(vector<Sequence> * db_seqs_to_update, vector<Sequence> * user_seqs_to_update);
    void retrieve_aligned_sequence_from_last_alignment_for_seq(Sequence * temp_seq);

    void load_sequences_from_last_alignment_into(vector<Sequence> & seqs);
    void add_seqs_from_db_to_seqs_vector(string alignname,vector<Sequence> * keep_seqs, vector<Sequence> & storedseqs);

public:
    SQLiteConstructor(string cn, vector <string> searchstr, bool searchlit, string genen, string genedb,
		      double mad_cut,double cover, double ident, string dbs,
		      string known_seq_filen, bool its, int numt,bool autom,
		      bool updb, string updf, bool ignoreleft, int shrinkthresh); // TODO: create a parameters struct that holds all this crap so we can just pass that
    ~SQLiteConstructor(){}
    void set_only_names_from_file(string filename, bool containshi, bool containswild);
    void set_exclude_names_from_file(string filename, bool, bool);
    void set_exclude_gi_from_file(string filename);
    void set_include_gi_from_file(string filename);
    void set_user_guide_tree(string filename, bool skipcheckdb);
    void set_user_fasta_file(string filename, bool skipcheckdb);
    void set_user_skip_search();
    void set_justseqquery(bool);
    int run(); //0 = success, 2 = update>20%
    string get_cladename();
    vector <string> get_search();
    string get_genename();
    string get_genedb();
    double get_madcutoff();
    double get_coverage();
    double get_identity();
    int get_numthreads();
    bool get_updatestatus();
    Tree * get_user_guide_tree_obj();
};


#endif	//_CONSTRUCTOR_H_
