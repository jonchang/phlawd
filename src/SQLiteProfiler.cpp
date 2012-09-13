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
 * SQLiteProfiler.cpp
 */

#include <string.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <unistd.h>

using namespace std;

#include "tree.h"
#include "node.h"
#include "tree_utils.h"
#include "tree_reader.h"
#include "sequence.h"
#include "fasta_util.h"

#include "utils.h"

#include "omp.h" 
#include "libsqlitewrapped.h"

#include "SQLiteProfiler.h"

#define UNDEF_DIST_SAME_ALN 666666
#define FAKE_ID 0
#define NOT_YET_PROFILED -1

template<class T>
inline std::string to_string(const T& t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

SQLiteProfiler::SQLiteProfiler(string gn, string gene_dbn, string cn,
		string dbs, bool autom, bool updb) :
		gene_name(gn), gene_db_fname(gene_dbn), clade_name(cn), source_db_fname(
				dbs), automated(autom), doing_update(updb) {
	use_orphan = false;
	using_guide_tree = false;
	user_guide_tree = NULL;
	final_alignment_dbid = -1;
	profile_dir_fname = gene_name + "_TEMPFILES/";
	gene_db = GeneDB(gene_db_fname);
}

void SQLiteProfiler::prelimalign() {
	mkdir(profile_dir_fname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH);
	bool standard = true;
	updated_original_alignment_dbids = vector<int>();
	if (doing_update == true) {
		vector<string> samefiles; //files were not updated at all
		bool newfiles = false; // if there are new files from MAD in align stage, kick out and tree as new for now
		vector<string> originalalnfiles;
		vector<string> curproffiles;
		gene_db.load_orig_alignment_names_into(originalalnfiles);
		gene_db.get_first_profile_alignments(curproffiles);
		if (originalalnfiles.size() != curproffiles.size()
				|| originalalnfiles.size() == 1)
			newfiles = true;
		else { //double check that there aren't new ones
			std::sort(originalalnfiles.begin(), originalalnfiles.end());
			std::sort(curproffiles.begin(), curproffiles.end());
			std::vector<string> v3;
			std::set_intersection(originalalnfiles.begin(),
					originalalnfiles.end(), curproffiles.begin(),
					curproffiles.end(), std::back_inserter(v3));
			if (v3.size() != originalalnfiles.size())
				newfiles = true;
		}
		if (newfiles == true) {
			gene_db.remove_profile_alignments();
			standard = true;
			doing_update = false;
		} else { //update the profiles
			standard = false;
			vector<string> updatedprofs;
			vector<int> updatedprofsnums; //these will be the profile ids that are updated
			vector<int> notupdatedprofsnums; //these will not be updated
			//TODO: start editing here, need the
			gene_db.get_updated_profs_names_delete_old(updatedprofs,
					updatedprofsnums, notupdatedprofsnums);
			updated_profile_alignment_dbids.clear();
			for (int i = 0; i < notupdatedprofsnums.size(); i++) {
				int tid = gene_db.get_deepest_profile_for_alignment(
						notupdatedprofsnums[i]);
				if (tid == -1) {
					updated_original_alignment_dbids.push_back(
							notupdatedprofsnums[i]);
					cout << "id: " << notupdatedprofsnums[i] << endl;
				} else if (count(updated_profile_alignment_dbids.begin(),
						updated_profile_alignment_dbids.end(), tid) == 0) {
					updated_profile_alignment_dbids.push_back(tid);
					cout << "tid: " << tid << endl;
				}
			}
		}
	}
	if (standard == true) { //not an update run

		// these are properties of the SQLiteProfiler class, accessed by various methods throughout
//		original_alignment_names = vector<string>();
//		original_alignment_dbids = vector<int>();

		cout << "getting alignments" << endl;
		gene_db.load_orig_alignment_names_into(original_alignment_names);

//		profile_id_name_map = map<int, string>();
		gene_db.migrate_original_alignments_and_load_new_profile_info_into(profile_id_name_map);
		gene_db.load_first_profile_ids_into(original_alignment_dbids);

	} else { // updaterun

		// only copy over the ones that are updated
		// only align those files that are updated
//		original_alignment_names = vector<string>();
//		original_alignment_dbids = vector<int>();
		cout << "getting updated alignment files" << endl;
//		profile_id_name_map = map<int, string>();

		gene_db.copy_alignments_to_first_profiles_updated(profile_id_name_map,
				updated_original_alignment_dbids);
		gene_db.load_orig_alignment_names_into(original_alignment_names);
		gene_db.load_first_profile_ids_into(original_alignment_dbids);
	}
}

void SQLiteProfiler::set_user_guide_tree(Tree * tree) {
	using_guide_tree = true;
	user_guide_tree = tree;
}

void SQLiteProfiler::run() {
	cout << "starting run" << endl;
	int finalaln;
	if (original_alignment_dbids.size() > 1) {
		map<string, string> numnames;
		map<string, string> namesnum;
		//	map<int, map<int,double> > original_alignment_distances;
		//TODO: make sure that this works with update , incomplete user guide tree
		if (using_guide_tree == true) {
			cout << "user guide tree" << endl;
			create_distances_user_tree(original_alignment_names, &numnames,
					&namesnum/*,&original_alignment_distances*/);
		} else {		//use ncbi tree
			cout << "ncbi guide tree" << endl;
			create_distances(clade_name/*,&original_alignment_distances*/);
		}
		//start profiling
		cout << "profiling" << endl;
		do_profile_alignments(/*original_alignment_distances*/);
//		finalaln
	} else {
		final_alignment_dbid = 1;
		clean_dbseqs(1);
	}

	if (doing_update == true)
		gene_db.toggle_updated_all_off();

	cout << "writing final alignment" << endl;
	rename_final_alignment(final_alignment_dbid);	//requires FINAL.aln
	//remove_outliers
	cout << "profile run completed" << endl;
}

void SQLiteProfiler::get_children(string in_id, vector<string> * in_ids,
		vector<string> * in_keepids) {
	Database conn(source_db_fname);
	string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id = " + in_id;
	sql += " and name_class = 'scientific name';";
	Query query(conn);
	query.get_result(sql);
//	StoreQueryResult R = query.store();
	while (query.fetch_row()) {
		string a;
		a = to_string(query.getstr());
		in_ids->push_back(a);
		in_keepids->push_back(a);
	}
}

vector<string> SQLiteProfiler::get_final_children(string id) {
	if (id == "1") {
		cout << "special case as id is root" << endl;
		Database conn(source_db_fname);
		vector<string> allids;
		string sql = "SELECT name,name_class,ncbi_id FROM taxonomy;";
		Query query(conn);
		query.get_result(sql);
		//StoreQueryResult R = query.store();
		while (query.fetch_row()) {
			//string tid = R[j][0].c_str();
			string tn = query.getstr();
			string cln = query.getstr();
			string ncbiid = query.getstr();
			if (cln.find("scientific") != string::npos
					&& tn.find("environmental") == string::npos
					&& cln.find("environmental") == string::npos) {
				allids.push_back(ncbiid); //was taxon id, now ncbi id
			}
		}
		query.free_result();
		return allids;
	} else {
		vector<string> ids;
		ids.push_back(id);
		vector<string> keepids;
		//check
		keepids.push_back(id);
		while (ids.size() > 0) {
			string tid = ids.back();
			ids.pop_back();
			get_children(tid, &ids, &keepids);
			//cout << ids.size() << endl;
		}
		vector<string> allids;
		vector<string> allnames;

		for (int i = 0; i < keepids.size(); i++) {
			Database conn(source_db_fname);
			string sql = "SELECT name FROM taxonomy WHERE ncbi_id = "
					+ keepids[i];
			sql += " and name_class = 'scientific name';";
			Query query(conn);
			query.get_result(sql);
			//StoreQueryResult R = query.store();
			while (query.fetch_row()) {
				string name;
				string cl;
				id = keepids[i];
				name = query.getstr();
//			cl = R[j][2].c_str();
//			size_t found1;
//			found1 = cl.find("environmental");
//			size_t found2;
//			found2 = cl.find("scientific");
//			if(found2!=string::npos && found1==string::npos){
				allids.push_back(id);
				allnames.push_back(name);
//			}
			}
		}
		allids.push_back(id);
		return allids;
	}
}

/*
 * replaces get_final_children and get_children with this
 * that uses the left and right values  to get all the children
 */

vector<string> SQLiteProfiler::get_left_right_children(string id) {
	vector<string> allids;
	Database conn(source_db_fname);
	string sql2 = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = "
			+ id;
	sql2 += " and name_class = 'scientific name';";
	Query query(conn);
	query.get_result(sql2);
	//StoreQueryResult R = query.store();
	string left, right;
	while (query.fetch_row()) {
		left = to_string(query.getval());
		right = to_string(query.getval());
	}

	string sql = "SELECT ncbi_id FROM taxonomy WHERE left_value >= " + left;
	sql += "AND right_value <= ";
	sql += right;
	sql += " and name_class = 'scientific name' ORDER BY left_value DESC;";
	Query query2(conn);
	query2.get_result(sql);
//	StoreQueryResult R2 = query2.store();
	while (query2.fetch_row()) {
		allids.push_back(to_string(query2.getval()));
	}
	return allids;
}

string SQLiteProfiler::get_right_one(vector<string> allids, Query & res) {
	while (res.fetch_row()) {
		int fd;
		string x;
		x = res.getstr();
		fd = count(allids.begin(), allids.end(), x);
		if (fd > 0)
			return x;
	}
	return "";
}

//ncbi one
void SQLiteProfiler::create_distances(
		string clade_name/*,map<int, map<int,double> > * original_alignment_distances*/) {
//    original_alignment_distances->clear();
	//get id for clade name
	// Make SQL string and execute it
	Database conn(source_db_fname);
//	StoreQueryResult R;
	Query query(conn);
	string sql;
	if (automated == false) {
		sql = "SELECT ncbi_id FROM taxonomy WHERE name = '" + clade_name
				+ "' and name_class = 'scientific name';";
		query.get_result(sql);
	} else if (automated == true) {
		sql = "SELECT ncbi_id FROM taxonomy WHERE ncbi_id = '" + clade_name
				+ "' and name_class = 'scientific name';";
		query.get_result(sql);
	}
	string cladeid;
	while (query.fetch_row()) {
		cladeid = to_string(query.getval());
	}
	cout << cladeid << endl;
	/*
	 * get left and right for the cladeid
	 */
	string sql2 = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = "
			+ cladeid;
	sql2 += " and name_class = 'scientific name';";
	Query query2(conn);
	query2.get_result(sql2);
//	StoreQueryResult R2 = query2.store();
	string cladeleft, claderight;
	while (query2.fetch_row()) {
		cladeleft = to_string(query2.getstr());
		claderight = to_string(query2.getstr());
	}
	/*
	 * end left and right
	 */
	vector<string> allids = get_final_children(cladeid);
	//vector<string> allids = get_left_right_children(cladeid);
	cout << "clade=" << cladeid << endl;
	//get the route to the clade name
//    for(int i=0;i<names.size();i++){
	map<int, string>::iterator it;
	for (it = profile_id_name_map.begin(); it != profile_id_name_map.end();
			it++) {
		Database conn(source_db_fname);
		//change to ncbi id
		sql = "SELECT ncbi_id FROM taxonomy WHERE ncbi_id = " + (*it).second
				+ ";";    //names[i]+";";
//	cout << sql << endl;
		Query query2(conn);
		query2.get_result(sql);
		string nameid = get_right_one(allids, query2);
		sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = " + nameid
				+ " and name_class = 'scientific name';";
		Query query4(conn);
		query4.get_result(sql);
		string parentid = get_right_one(allids, query4);
		vector<string> route;
		while (parentid != cladeid) {
			route.push_back(parentid);
			nameid = parentid;
			sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "
					+ nameid + " and name_class = 'scientific name';";
			//cout << sql << endl;
			Query query5(conn);
			query5.get_result(sql);
			parentid = get_right_one(allids, query5);
		}
		route.push_back(parentid);

		map<int, double> tdistance;
		map<int, string>::iterator it2;
//	for(int j=0;j<names.size();j++){
		for (it2 = profile_id_name_map.begin();
				it2 != profile_id_name_map.end(); it2++) {
			if (it2 != it) {
				//sql = "SELECT ncbi_id FROM taxonomy WHERE name = '"+names[j]+"'  and name_class = 'scientific name';";
				//change to ncbi id
				sql = "SELECT ncbi_id FROM taxonomy WHERE ncbi_id = "
						+ (*it2).second + ";";                //names[j]+";";
				Query query5(conn);
				query5.get_result(sql);
				string jnameid = get_right_one(allids, query5);
				double distance = 0;
				while ((int) count(route.begin(), route.end(), jnameid) == 0) {
					distance += 1;
					sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "
							+ jnameid + " and name_class = 'scientific name';";
					Query query6(conn);
					query6.get_result(sql);
					jnameid = get_right_one(allids, query6);
				}
				for (int k = 0; k < route.size(); k++) {
					if (route[k] == jnameid)
						distance += k;
				}
				tdistance[(*it2).first] = distance;
			} else {
				tdistance[(*it2).first] = UNDEF_DIST_SAME_ALN;
			}
		}
		(/***/original_alignment_distances)[(*it).first] = tdistance;
		cout << "distances complete: " << (*it).second << " " << (*it).first
				<< endl;
	}
}

//user tree one
void SQLiteProfiler::create_distances_user_tree(vector<string> ifile_names,
		map<string, string> * numnames,
		map<string, string> * namesnum/*, map<int, map<int,double> > * original_alignment_distances*/) {
	//get the list of nodes for which distances are required
	vector<Node *> nodesfordist;
	for (int i = 0; i < ifile_names.size(); i++) {
		for (int j = 0; j < user_guide_tree->getNodeCount(); j++) {
			if (user_guide_tree->getNode(j)->getName() == ifile_names[i])
				nodesfordist.push_back(user_guide_tree->getNode(j));
		}
	}
	for (int i = 0; i < nodesfordist.size(); i++) {
		vector<double> tdistance;
		for (int j = 0; j < nodesfordist.size(); j++) {
			if (i == j) {
				tdistance.push_back(UNDEF_DIST_SAME_ALN);
			} else {
				double distance = get_distance_between_two_nodes(
						user_guide_tree, nodesfordist[i], nodesfordist[j]);
				tdistance.push_back(distance);
			}
		}
		std::ostringstream stm;
		stm << i;
		numnames->insert(
				pair<string, string>(stm.str(), nodesfordist[i]->getName()));
		namesnum->insert(
				pair<string, string>(nodesfordist[i]->getName(), stm.str()));
//	original_alignment_distances->push_back(tdistance);
	}
}

void SQLiteProfiler::find_next_original_alignment_set_to_profile(
		int * best_seed_alignment, vector<int> * closest_matches) {

	/* This function uses the alignments from original_alignments_to_profile as
	 * seeds to identify the taxonomically closest set of *all* original alignments
	 * that includes at least one alignment from original_alignments_to_profile. It
	 * stores the db id of the first such alignment from original_alignments_to_profile
	 * in best_seed_alignment, and the set of its closest matches (which may include
	 * alignments that have already been removed from original_alignments_to_profile)
	 * in closest_matches.
	 */

	map<int, double> best_seed_distances;
	double shortest_distance = 10000;

	for (int i = 0; i < original_alignments_to_profile.size(); i++) {
		int aln1 = original_alignments_to_profile.at(i);

		for (int j = 0; j < original_alignment_dbids.size(); j++) {
			int aln2 = original_alignment_dbids.at(j);

			// compare every to-be-profiled original aln to all original alns
			double d = original_alignment_distances[aln1][aln2];
			if (d < shortest_distance && d != UNDEF_DIST_SAME_ALN) {

				// remember the to-be-profiled alignment with smallest distance to some other aln
				shortest_distance = d;
				*best_seed_alignment = aln1;
			}
		}
	}

	best_seed_distances = original_alignment_distances[*best_seed_alignment];
	closest_matches->clear();

	// find all the alignments that are equally closest to *best_seed_alignment
	for (int j = 0; j < original_alignment_dbids.size(); j++) {
		int alt_sib = original_alignment_dbids.at(j);
		double d = best_seed_distances[alt_sib];
		if (d == shortest_distance && d != UNDEF_DIST_SAME_ALN)
			closest_matches->push_back(alt_sib);
	}
}

void SQLiteProfiler::clean_alignment(string filename) {

	/* just cleaning the file that is there
	 * 0.1 is the current limit
	 * percent should be, if you are missing more than this, then remove, so 0.9 is now and 0.5 is genome removal
	 */

	double threshold = 0.9;

	// read the alignment into tempalseqs
	FastaUtil fu;
	vector<Sequence> tempalseqs;
	fu.readFile(profile_dir_fname + filename, tempalseqs);

	cout << "cleaning seqs for " << filename << endl;
	int seqlength = tempalseqs[0].get_sequence().size();
	float nseqs = float(tempalseqs.size());
	vector<int> cols_to_remove;

	// for each column in the alignment
	for (int j = 0; j < seqlength; j++) {
		int gaps = 0;
		char cell;
		// get the number of sequences with missing data
		for (int i = 0; i < tempalseqs.size(); i++) {
			cell = tempalseqs[i].get_sequence()[j];
			if (cell == '-' || cell == 'N' || cell == 'n')
				gaps += 1;
		}
		// if we exceed the threshold, mark this column for removal
		double prop_missing = gaps / nseqs;
		if (prop_missing > threshold)
			cols_to_remove.push_back(j);
	}

	// for each sequence
	for (int i = 0; i < tempalseqs.size(); i++) {
		string a;
		// for each column
		for (int j = 0; j < seqlength; j++) {
			// record all columns not in the exclusion list
			if (count(cols_to_remove.begin(), cols_to_remove.end(), j) == 0)
				a += tempalseqs[i].get_sequence()[j];
		}
		// record cleaned sequence
		tempalseqs[i].set_sequence(a);
	}

	// remove the input (dirty) file, replace it with the clean one
	remove((profile_dir_fname + filename).c_str());
	fu.writeFileFromVector(profile_dir_fname + filename, tempalseqs);
}

void SQLiteProfiler::clean_dbseqs(int alignid) {
	double percent = 0.9;
	vector<Sequence> tempalseqs;
	gene_db.get_profile_align_seqs(alignid, tempalseqs);
	cout << "cleaning seqs" << endl;
	int seqlength = tempalseqs[0].get_aligned_seq().size();
	float fseql = float(tempalseqs.size());
	vector<int> removeem;
	for (int j = 0; j < seqlength; j++) {
		int gaps = 0;
		for (int i = 0; i < tempalseqs.size(); i++) {
			if (tempalseqs[i].get_aligned_seq()[j] == '-'
					|| tempalseqs[i].get_aligned_seq()[j] == 'N'
					|| tempalseqs[i].get_aligned_seq()[j] == 'n')
				gaps += 1;
		}
		double curp = gaps / fseql;
		if (curp > percent) {
			removeem.push_back(j);
		}
	}
	for (int i = 0; i < tempalseqs.size(); i++) {
		string a;
		for (int j = 0; j < seqlength; j++) {
			if (count(removeem.begin(), removeem.end(), j) == 0)
				a += tempalseqs[i].get_aligned_seq()[j];
		}
		tempalseqs[i].set_aligned_seq(a);
	}
	gene_db.update_profile_align_seqs(alignid, tempalseqs);
}

/*
 * TODO: fix the orphans
 */
void SQLiteProfiler::do_profile_alignments(/*map<int, map<int,double> > original_alignment_distances*/) {

	/* This does all the profile alignments.
	 *
	 * We keep track of original and profile alignments in two vectors:
	 * original_alignments_to_profile and intermediate_profile_alignments,
	 * which store the database ids of the alignments. Note that the
	 * original alignments in the database itself have already been copied
	 * into the profile_alignments table, so all the dbids in both of these
	 * vectors represent rows in the profile_alignments table in the db.
	 *
	 * We start with the closest pair of original alignments. From these we
	 * create a profile alignment, store it in intermediate_profile_alignments,
	 * and remove the original alignments from original_alignments_to_profile.
	 * Then we find the closest pair of remaining original alignments, and
	 * repeat this process until no original alignments remain. After this we
	 * cross-align all the profile alignments, adding each newly created profile
	 * to intermediate_profile_alignments and removing both its children,
	 * until only one profile alignment remains, which is our final alignment.
	 */

	bool muscle = true;

	// commented this out because the log file never seems to be used
	//  cout << "writing everything to record.log" << endl;
	//	string recordname = profile_dir_fname + "/record.log";
	//	ofstream ofs(recordname.c_str());

	if (doing_update) {

		original_alignments_to_profile = updated_original_alignment_dbids;
		intermediate_profile_alignments = updated_profile_alignment_dbids;

		cout << "original alignments remaining to profile: " << endl;
		for (int i = 0; i < original_alignments_to_profile.size(); i++)
			cout << original_alignments_to_profile[i] << " ";
		cout << endl;

		cout << "intermediate profile alignments remaining to cross-align: ";
		for (int i = 0; i < intermediate_profile_alignments.size(); i++)
			cout << intermediate_profile_alignments[i] << " ";
		cout << endl;

	} else
		// not doing update; grab all the original alignments
		original_alignments_to_profile.assign(original_alignment_dbids.begin(),
				original_alignment_dbids.end());

	int current_matched_alignment;// reused by each iteration of the upcoming while loop
	int last_profile_completed; // remember the last profile in case it is the final one

	while (original_alignments_to_profile.size() > 0) { // first, the original alignments

		// find the closest equidistant set of remaining original alignments
		vector<int> * all_closest_matches = new vector<int>();
		find_next_original_alignment_set_to_profile(&current_matched_alignment,
				all_closest_matches);
		int n_closest_matches = all_closest_matches->size();

		// DEBUG
		cout << "next original alignment to profile: "
				<< current_matched_alignment << endl;
		for (int i = 0; i < all_closest_matches->size(); i++) {
			cout << "\tmatched " << all_closest_matches->at(i) << endl;
		}
		cout << "removing " << current_matched_alignment
				<< " from original alignments to profile" << endl;

		// don't profile current_matched_alignment against itself
		remove_from_original_alignments_to_profile(current_matched_alignment);

		if (n_closest_matches == 1) {
			int best_original_match = all_closest_matches->at(0);

			// will be the actual alignment that we use in the new profile
			int best_match;
			bool best_match_is_profiled;

			int best_profiled_match = gene_db.get_deepest_profile_for_alignment(
					best_original_match);
			if (best_profiled_match == NOT_YET_PROFILED) { // use the original alignment itself
				best_match = best_original_match;
				best_match_is_profiled = false;

			} else { // otherwise use the profile
				best_match = best_profiled_match;
				best_match_is_profiled = true;
			}

			// make the new profile
			cout << "profiling " << current_matched_alignment << " and "
					<< best_match << endl;
			int p = make_new_profile_alignment(current_matched_alignment,
					best_match);

			// remove best_match so we don't profile it again
			if (best_match_is_profiled) {
				cout << "removing " << best_match
						<< " from intermediate profiles" << endl;
				remove_from_intermediate_profiles(best_match);
			} else {
				cout << "removing " << best_match
						<< " from original alignments to profile" << endl;
				remove_from_original_alignments_to_profile(best_match);
			}

			last_profile_completed = p;

		} else { // there is more than one equally closest alignment to current_matched_alignment

			int first_profile_match = NOT_YET_PROFILED;
//            bool found_profiled_match = false;

			// walk the alignments in the list of equally closest matches...
			for (int j = 0; j < n_closest_matches; j++) {
				// ...checking each to see if it has already been profiled
				int m = all_closest_matches->at(j);
				int this_match = gene_db.get_deepest_profile_for_alignment(m);
				if (this_match != NOT_YET_PROFILED) {
					// if we find one, stop looking
					first_profile_match = this_match;
//                    found_profiled_match = true;
					break;
				}
			}

			if (first_profile_match != NOT_YET_PROFILED) { // we found a profile to align to current_matched_alignment

				// do the profile alignment, and remove the intermediate profile child alignment
				cout << "profiling " << current_matched_alignment << " and "
						<< first_profile_match << endl;
				int p = make_new_profile_alignment(current_matched_alignment,
						first_profile_match);
				cout << "removing " << first_profile_match
						<< " from intermediate profiles" << endl;
				remove_from_intermediate_profiles(first_profile_match);
				last_profile_completed = p;

			} else { // we did not find a profiled alignment in the set of closest matches

				int best_original_match;
				double best_score = -1;

				// find the alignment that scores highest against current_matched_profile
				for (int j = 0; j < n_closest_matches; j++) {
					int this_match = all_closest_matches->at(j);
					make_muscle_profile(current_matched_alignment, this_match);
					double score = get_muscle_spscore("TEMPOUT.PROFILE");
					if (score > best_score) {
						best_score = score;
						best_original_match = this_match;
					}
				}
//                current_closest_match = best_match;

				// do the profile alignment, and remove the original child alignment
				cout << "profiling " << current_matched_alignment << " and "
						<< best_original_match << endl;
				int p = make_new_profile_alignment(current_matched_alignment,
						best_original_match);
				cout << "removing " << best_original_match
						<< " from original alignments to profile" << endl;
				remove_from_original_alignments_to_profile(best_original_match);
				last_profile_completed = p;
			}
		}
		delete all_closest_matches;
	}

	cout << "processing profile files" << endl;
	for (int i = 0; i < intermediate_profile_alignments.size(); i++) {
		cout << intermediate_profile_alignments[i] << endl;
	}

	while (intermediate_profile_alignments.size() > 1) {
		// TODO: choose the best, for now just choose any two
		int aln1 = intermediate_profile_alignments[0];
		int aln2 = intermediate_profile_alignments[1];

		int p = make_new_profile_alignment(aln1, aln2);
		remove_from_intermediate_profiles(aln1);
		remove_from_intermediate_profiles(aln2);

		last_profile_completed = p;
	}

	final_alignment_dbid = last_profile_completed;
}

int SQLiteProfiler::make_new_profile_alignment(int aln1, int aln2) {

	// create new db record for the upcoming profile alignment
	int id = gene_db.add_empty_intermediate_profile(aln1, aln2);

	// record the new profile id
	intermediate_profile_alignments.push_back(id);

	// do the alignment and store it in the db
	make_muscle_profile(aln1, aln2);
	clean_alignment("TEMPOUT.PROFILE");
	match_and_add_profile_alignment_to_db(id);

	return id;
}

void SQLiteProfiler::remove_from_intermediate_profiles(int value) {
	remove_int_from_vector(value, intermediate_profile_alignments);
}

void SQLiteProfiler::remove_from_original_alignments_to_profile(int value) {
	remove_int_from_vector(value, original_alignments_to_profile);
}

void SQLiteProfiler::remove_int_from_vector(int value, vector<int> & v) {

	vector<int>::iterator index;
	index = find(v.begin(), v.end(), value);

	if (index == v.end())
		cerr << value << " not found" << endl;

	v.erase(index);

}

void SQLiteProfiler::make_muscle_profile(int aln1, int aln2) {

	/* accepts the db ids for two preexisiting alignments (in the
	 * profile_alignments table), writes each to a file, and feeds
	 * these files to muscle to profile align them. the resulting
	 * alignment is stored in TEMPOUT.PROFILE. */

	// delete any previous temp alignment files
	remove((profile_dir_fname + "TEMP1.PROFILE").c_str());
	remove((profile_dir_fname + "TEMP2.PROFILE").c_str());
	remove((profile_dir_fname + "TEMPOUT.PROFILE").c_str());

	// write the child alignments to temp files
	gene_db.write_profile_alignment_to_file(aln1,
			profile_dir_fname + "TEMP1.profile");
	gene_db.write_profile_alignment_to_file(aln2,
			profile_dir_fname + "TEMP2.profile");

	// build the muscle call
	string cmd = "muscle -profile -maxmb 5000 -in1 ";
	cmd += profile_dir_fname + "TEMP1.profile -in2 ";
	cmd += profile_dir_fname + "TEMP2.profile -out ";
	cmd += profile_dir_fname + "TEMPOUT.PROFILE 2> ";
	cmd += profile_dir_fname + "muscle.out";

	cout << "aligning" << endl;
	cout << cmd << endl;

	// call muscle to do profile alignment
	system(cmd.c_str());

	// if new profile alignment was not created, test will fail and program will exit
	validate_outfile_exists(profile_dir_fname + "TEMPOUT.PROFILE");

}

string SQLiteProfiler::get_name_from_tax_id(string taxid) {
	Database conn(source_db_fname);
	string sql = "SELECT name FROM taxonomy WHERE ncbi_id = " + taxid
			+ " AND name_class = 'scientific name'";
	Query query(conn);
	query.get_result(sql);
//	StoreQueryResult R = query.store();
	string ret;
	while (query.fetch_row()) {
		string a;
		a = query.getstr();
		ret = a;
	}
	return ret;
}

/*
 * rename the FINAL.aln.cln file to FINAL.aln.cln.rn using the table in ITS.gi
 */
void SQLiteProfiler::rename_final_alignment(int alignid) {
	string fn1 = gene_name + ".FINAL.aln.rn";
	gene_db.write_profile_alignment_with_names_to_file(alignid, fn1, false);
	fn1 = gene_name + ".FINAL.aln";
	gene_db.write_profile_alignment_with_names_to_file(alignid, fn1, true);
}

double SQLiteProfiler::get_muscle_spscore(string filename) {
	remove("prolog");
	string cmd = "muscle -maxmb 5000 -spscore ";
	cmd += profile_dir_fname;
	cmd += filename;
	cmd += " -log ";
	cmd += profile_dir_fname + "prlog 2> ";
	cmd += profile_dir_fname + "muscle.out";
	cout << cmd << endl;
	/*    FILE *fp = popen(cmd.c_str(), "r" );
	 char buff[1000];
	 while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	 string line(buff);
	 }
	 pclose( fp );
	 */
	system(cmd.c_str());
	//read file
	double score = 0;
	ifstream ifs((profile_dir_fname + "prlog").c_str());
	string line;
	while (getline(ifs, line)) {
		TrimSpaces(line);
		size_t se = line.find("="); // Find the first character position from reverse af
		// if all spaces or empty return an empty string
		if (se != string::npos) {
			vector<string> tokens;
			string del("=");
			Tokenize(line, tokens, del);
			score = atof(tokens[2].c_str());
		}
	}
	ifs.close();
	return score;
}

void SQLiteProfiler::validate_outfile_exists(string filename) {
	/* test whether filename exists, if not, the
	 * program will exit with an error */

	struct stat filestatus;
	stat(filename.c_str(), &filestatus);
	cout << filestatus.st_size << " bytes\n";
	if (int(filestatus.st_size) <= 64) {
		cerr << "problem: empty file " << filename << " created" << endl;
		exit(0);
	}
}

/*
 * assumes that the id of the seq is the sqlite id 
 */
void SQLiteProfiler::match_and_add_profile_alignment_to_db(int profileid) {
	FastaUtil fu;
	vector<Sequence> tempalseqs;
	fu.readFile(profile_dir_fname + "TEMPOUT.PROFILE", tempalseqs);
	vector<Sequence> seqs;
	for (int i = 0; i < tempalseqs.size(); i++) {
		Sequence tseq(tempalseqs[i].get_id(), tempalseqs[i].get_sequence());
		seqs.push_back(tseq);
	}
	cout << "adding sequences to profile alignment table" << endl;
	gene_db.add_sequences_for_profile_alignment(profileid, seqs);
}

//THIS WAS THE OLD quicktree based outlier removal
/*
 void SQLiteProfiler::calculate_for_removal_quicktree(vector<Sequence> * seqs,
 map<string,double> & allmeans){
 FastaUtil seqwriter;
 const string fn1 = "TEMPFILES/tempremoval";
 seqwriter.writeFileFromVector(fn1,*seqs);

 string phcmd = "phyutility -concat -in TEMPFILES/tempremoval -out TEMPFILES/outfile.nex";
 FILE *phfp = popen(phcmd.c_str(), "r" );
 pclose( phfp );

 cout << phcmd << endl;

 ifstream infile;
 ofstream outfile;
 infile.open ("TEMPFILES/outfile.nex",ios::in);
 outfile.open ("TEMPFILES/outfile.stoc",ios::out);
 bool begin = false;
 bool end = false;
 string line;
 // convert to stockholm format
 while(getline(infile,line)){
 if (line.find("MATRIX") != string::npos){
 begin = true;
 }else if ((begin == true && end == false) && line.find_first_of(";") != string::npos){
 end = true;
 }else if (begin == true && end == false){
 std::string::size_type begin = line.find_first_not_of("\t");
 //std::string::size_type end   = line.find_last_not_of("\t");
 std::string::size_type end = line.size();
 std::string trimmed = line.substr(begin, end-begin + 1);
 outfile << trimmed << endl;
 }
 }
 infile.close();
 outfile.close();

 const char * cmd = "quicktree -in a -out m TEMPFILES/outfile.stoc > TEMPFILES/dist";
 cout << "calculating distance" << endl;
 FILE *fp = popen(cmd, "r" );
 char buff[1000];
 while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
 string line(buff);
 }
 pclose( fp );

 //new
 vector<string> ids;
 vector<vector<double> > nums;
 bool first = true;

 ifstream pfile ("TEMPFILES/dist");
 vector<string> tokens;
 int nspecies;
 if (pfile.is_open()){
 int curspecies = 0;
 while (! pfile.eof() ){
 // record the id as the first token
 // and the numbers below as the numbers
 if (first == false){
 getline (pfile,line);
 string del(" \t");
 tokens.clear();
 Tokenize(line, tokens, del);
 if(tokens.size() >= 1){
 ids[curspecies] = tokens.at(0);
 double n1;
 for(int j = curspecies; j < nspecies;j++){
 n1 = atof(tokens.at(j+1).c_str());
 nums[curspecies][j] = n1;
 nums[j][curspecies] = n1;
 }
 curspecies += 1;
 }
 }else{
 first = false;
 getline (pfile,line);
 TrimSpaces(line);
 nspecies = atoi(line.c_str());
 vector<double> cols(nspecies, 0);
 nums = vector< vector<double> >(nspecies, cols);
 ids = vector<string>(nspecies);
 }
 }
 pfile.close();
 }
 // calculate the means
 vector<double> mns(ids.size());
 for(int i=0;i<nums.size();i++){
 allmeans[ids[i]] = mean(nums[i]);
 }
 }

 void SQLiteProfiler::remove_outliers(){
 FastaUtil seqreader;
 vector<Sequence> * sequences = new vector<Sequence>();
 seqreader.readFile(profile_dir_fname+"FINAL.aln", *sequences);
 int numseqs = sequences->size();
 cout << numseqs << endl;
 map<string,double> allmeans;
 if(numseqs < 5000){
 calculate_for_removal_quicktree(sequences,allmeans);
 }else{
 int NBREAKS = 10;
 for (int i=0;i<NBREAKS;i++){
 vector<Sequence> tempsc1;
 if((i+1) < NBREAKS){
 for(int j=(i*(numseqs/NBREAKS));j < ((numseqs/NBREAKS)*(i+1));j++){
 //TODO : there was a pointer problem here
 Tempsc1.push_back(sequences->at(j));
 }
 }else{
 for(int j=(i*(numseqs/NBREAKS));j < numseqs;j++){
 //TODO : there was a pointer problem here
 tempsc1.push_back(sequences->at(j));
 }
 }
 calculate_for_removal_quicktree(&tempsc1,allmeans);
 }
 }
 // calculate the means
 vector<double> mns(numseqs);
 for(int i=0;i<numseqs;i++){
 //TODO : there was a pointer problem here
 mns[i] = allmeans[sequences->at(i).get_id()];
 }
 double sd = stdev(mns);
 double mn = mean(mns);
 double dev = sd*1.5+mn;//obviously changeable
 // read in the fasta file
 vector<Sequence> sc1;
 for(int i=0;i<mns.size();i++){
 if (mns[i] < dev){
 sc1.push_back(sequences->at(i));
 }else{
 //cout << i << endl;
 }
 }
 FastaUtil seqwriter;
 const string fn1 = profile_dir_fname+"FINAL.aln.cln";
 seqwriter.writeFileFromVector(fn1,sc1);
 }
 */
