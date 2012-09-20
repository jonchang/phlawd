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
 * SQLiteProfiler.h
 */

#ifndef SQLITEPROFILER_H_
#define SQLITEPROFILER_H_

#include "tree.h"
#include "genedb.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

#include "libsqlitewrapped.h"

class SQLiteProfiler {
private:

	bool use_orphan;								// ?
	bool automated;									// ?
	bool doing_update;								// are we performing an update
	bool using_guide_tree;							// do we have a user-supplied guide tree

	int final_alignment_dbid;						// the fully-inclusive profile alignment

	string gene_name;								// name of locus that phlawd is aligning (from config file)
	string clade_name;								// search clade (from config file)
	string gene_db_fname;							// database to hold phlawd results
	string source_db_fname;							// source database
	string profile_dir_fname;						// location of profile alignments

	Tree * user_guide_tree;							// the user supplied guide tree, if it exists
	GeneDB gene_db;									// the interface object for the phlawd results db

	vector<string> original_alignment_names;		// ?
	vector<int> original_alignment_dbids;			// db ids of the original alignments (not profiles)

	vector<int> original_alignment_dbids_to_update;		// used if we are doing an update run
	vector<int> profile_alignment_dbids_to_update;		// used if we are doing an update run

	vector<int> intermediate_profile_alignments;	// db ids of the profile alignments remaining to be cross-profiled
	vector<int> original_alignments_to_profile;		// db ids of the original alignments remaining to be profiled

	map<int, string> profile_id_name_map;			// hashmap associating profile alignments with their database ids

	// taxonomic distances between every pair of original alignments,
	// indexed by each alignment's db id in the profile_alignments table
	map<int, map<int,double> > original_alignment_distances;

	//functions
	void get_children(string in_id, vector<string> * in_ids, vector<string> * in_keepids);
	vector<string> get_final_children(string id);
	string get_right_one(vector<string> allids, Query & res);
	vector<string> get_left_right_children(string id);
	void create_distances(/*string clade_name*/);
	void create_distances_user_tree(vector<string> names, map<string, string> * numnames, map<string, string> * namesnum);
	void find_next_original_alignment_set_to_profile(int * matched_alignment, vector<int> * closest_matches);
	void clean_aligned_sequences(vector<Sequence> & seqs, float site_threshold, float seq_threshold);
	void clean_alignment(string infile);
	void clean_dbseqs(int alignid);
	void do_profile_alignments();
	string get_name_from_tax_id(string taxid);
	void write_final_alignments(int alignid);
	void make_muscle_profile(int aln1, int aln2);
	double get_muscle_spscore(string filename);

	void validate_file(string & filename); // TODO: this should probably go in some other class, like util or something

	void match_and_add_profile_alignment_to_db(int profileid);

	void remove_from_intermediate_profiles(int value);
	void remove_from_original_alignments_to_profile(int value);
	void remove_int_from_vector(int value, vector<int> & the_vector);

	int make_new_profile_alignment(int aln1, int aln2);

public:
	SQLiteProfiler(string gn, string gene_dbn, string cn, string dbs,
			bool autom, bool updb);
	void prelimalign();
	void run();
	void set_user_guide_tree(Tree * tree);
};

#endif /* MQPROFILER_H_ */
