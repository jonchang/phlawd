#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <stack>
#include <algorithm>

#include "genedb.h"
#include "libsqlitewrapped.h"
#include "utils.h"
#include "fasta_util.h"
#include "sequence.h"

#define NOT_YET_PROFILED -1;

using namespace std;

template<class T>
inline std::string to_string(const T& t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

GeneDB::GeneDB() {
}

GeneDB::GeneDB(string nm) :
		name(nm) {

//	cout << *db << endl;
//	cout << &db << endl;
//	cout << *db << endl;
//	cout << &db << endl;

	/*	db = &*temp_handle;

	 if (rc != SQLITE_OK) {
	 cerr << "Could not open sqlite3 connection to '" << name.c_str() <<  "'; error " << rc << endl;
	 sqlite3_close(db);
	 } */

}

GeneDB::~GeneDB() {
//	sqlite3_close(db);
}
/*
 * Should check to see if the database exists, if not 
 * then it will initialize it. If so then it will make
 * sure that the structure is correct.
 */
/*
 * SQL structure as of 0.3a
 * 
 * alignments
 * -
 * 
 * sequences
 * -
 *
 * sequence_alignment_map
 * -
 * 
 *
 * alignment_alignment_map
 * -
 * 
 * record
 * -
 */

/******************************************************************************************
 *
 *	the following private functions define interaction with the database itself, as
 *	much as that is possible.
 *
 *	haven't yet found a way to open a single database connection that can be shared by
 *	different member functions, so we still do that individually within each public
 *	function that calls these.
 *
 */

string GeneDB::get_name() {

	/* just return the name of the phlawd db we are working with */
	return name;
}

void GeneDB::execute_query(sqlite3 * db, const char * query) {

	/* this function actually calls the sqlite3 library functions,
	 * which only accept a c-style string for the query. it is
	 * called by an identically-named function that is provided
	 * for convenience, which accepts queries as c++ strings, and
	 * just converts them to c-style before passing them here.
	 */

	int rc;
	char *zErrMsg = 0;

	rc = sqlite3_exec(db, query, 0, 0, &zErrMsg);
	if (rc != SQLITE_OK) {
		cerr << "Could not execute query: " << sqlite3_errmsg(db) << endl << "The query was: " << query << endl;
		abort();
	}
}

void GeneDB::execute_query(sqlite3 * db, string & query) {

	/* this function is just a wrapper for the c-style string
	 * implementation of execute_query, so that we can provide
	 * either c or c++ formatted strings without worry.
	 */

	const char * query_c_str = query.c_str();
	execute_query(db, query_c_str);
}

void GeneDB::begin_transaction(sqlite3 * db) {
	execute_query(db, "BEGIN TRANSACTION");
}

void GeneDB::commit_transaction(sqlite3 * db) {
	execute_query(db, "COMMIT TRANSACTION");
}

void GeneDB::execute_simple_transaction(sqlite3 * db, string & sql) {
	begin_transaction(db);
	execute_query(db, sql);
	commit_transaction(db);
}

void GeneDB::initialize(bool overwrite) {

	/* deletes the existing database (if it exists, and if
	 * overwriting is allowed), creates a new database, and
	 * installs the requisite database structure.
	 */

	bool install_new = true; 		// guilty until proven innocent

	ifstream ifile(name.c_str()); 	// check if db file already exists
	if (ifile) {
		if (overwrite) {
			// erase the existing db if we're allowed
			remove(name.c_str());
			install_new = true;
		} else { // bail
			cerr << name << " exists" << endl;
			install_new = false;

			// TODO: if the db exists, and we aren't going to erase it, we
			// should to make sure that it has the right structure
		}
	}

	if (install_new) {
		sqlite3 * db;
		if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
			cout << "connected to " << name << endl;

			begin_transaction(db);
			execute_query(db, "create table alignments(id INTEGER PRIMARY KEY, alignname VARCHAR(255),updated INTEGER);");
			execute_query(db, "create index alignments_alignname on alignments(alignname);");
			execute_query(db, "create table profile_alignments(id INTEGER PRIMARY KEY, child1 INTEGER, child2 INTEGER, name VARCHAR(255));");
			// make index
			execute_query(db, "create table sequence_profile_map(id INTEGER PRIMARY KEY, sequence_id INTEGER, profile_id INTEGER, sequence TEXT);");
			// make index
			execute_query(db, "create table sequences(id INTEGER PRIMARY KEY, ncbi_id INTEGER, accession VARCHAR(255), name VARCHAR(255), sequence TEXT);");
			execute_query(db, "create index sequences_ncbi_id on sequences(ncbi_id);");
			execute_query(db, "create index sequences_name on sequences(name);");
			execute_query(db, "create table sequence_alignment_map(id INTEGER PRIMARY KEY, sequence_id INTEGER, alignment_id INTEGER, sequence TEXT);");
			execute_query(db, "create index sequence_alignment_map_sequence_id on sequence_alignment_map(sequence_id);");
			execute_query(db, "create index sequence_alignment_map_alignment_id on sequence_alignment_map(alignment_id);");
			commit_transaction(db);

			cout << "database created" << endl;

		} else {
			cout << "could not create the database: " << name << endl;
		}
		sqlite3_close(db);
	}
}

int GeneDB::insert_sequence(sqlite3 * db, Sequence & seq) {

	/* inserts a row in the db 'sequences' table containing
	 * the information from a Sequence object, and returns the
	 * corresponding db record id for the new row.
	 */

	string sql = "insert into sequences (ncbi_id, accession, name, sequence) values (";
	sql += seq.get_ncbi_tax_id() + ",'";
	sql += seq.get_ncbi_gi_id() + "','";
	sql += seq.get_name() + "','";
	sql += seq.get_sequence() + "');";

	execute_query(db, sql);

	return sqlite3_last_insert_rowid(db);
}

int GeneDB::insert_alignment(sqlite3 * db, string & alignment_name, int updated) {

	/* inserts a row in the db 'alignments' table containing
	 * the information provided as arguments, and returns the
	 * corresponding db record id for the new row.
	 */

	string sql = "insert into alignments (alignname, updated) values ('";
	sql += alignment_name + "', ";
	sql += to_string(updated) + ");";
	execute_query(db, sql);

	return sqlite3_last_insert_rowid(db);
}

int GeneDB::insert_profile_alignment(sqlite3 * db, string & profile_name, int child1_id, int child2_id) {

	/* this is the private method used by other methods within GeneDB.
	 * it inserts a row in the db 'profile_alignments' table containing
	 * the information provided as arguments, and returns the
	 * corresponding db record id for the new row.
	 */

	string sql = "insert into profile_alignments (name, child1, child2) values ('";
	sql += profile_name + "',";
	sql += to_string(child1_id) + ",";
	sql += to_string(child2_id) + ");";
	execute_query(db, sql);

	return sqlite3_last_insert_rowid(db);
}

int GeneDB::insert_sequence_alignment_map(sqlite3 * db, int sequence_id, int alignment_id, string & sequence) {

	/* inserts a row in the db 'sequence_alignment_map' table
	 * containing the information provided as arguments, and
	 * returns the corresponding db record id for the new row.
	 */

	string sql = "insert into sequence_alignment_map (sequence_id, alignment_id, sequence) values (";
	sql += to_string(sequence_id) + ",";
	sql += to_string(alignment_id) + ",'";
	sql += sequence + "');";
	execute_query(db, sql);

	return sqlite3_last_insert_rowid(db);
}

int GeneDB::insert_sequence_profile_map(sqlite3 * db, int sequence_id, int profile_id, string & sequence) {

	/* inserts a row in the db 'sequence_profile_map' table
	 * containing the information provided as arguments, and
	 * returns the corresponding db record id for the new row.
	 */

	string sql = "insert into sequence_profile_map (sequence_id, profile_id, sequence) values (";
	sql += to_string(sequence_id) + ",";
	sql += to_string(profile_id) + ",'";
	sql += sequence + "');";
	execute_query(db, sql);

	return sqlite3_last_insert_rowid(db);
}

void GeneDB::add_sequence_to_alignment(int alignment_id, Sequence & sequence_obj) {

	/* calls insert_sequence_alignment_map to create a new
	 * record associating the alignment id and the sequence
	 * contained in the object passed as parameters.
	 *
	 * assumes that the passed sequence object has already
	 * been stored in the db and contains a valid database
	 * id. if not, this could do bad things.
	 */

	// gather the info to be added to the db
	int sequence_id = sequence_obj.get_sqlite_id();
	string sequence = sequence_obj.get_aligned_seq();

	// do the insert, or fail if we can't open the db
	sqlite3 *db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		begin_transaction(db);
		insert_sequence_alignment_map(db, sequence_id, alignment_id, sequence);
		commit_transaction(db);
	} else
		cerr << "Could not add sequence " << sequence_id << " because the database " << name << " could not be opened";

	sqlite3_close(db);

}

void GeneDB::add_sequences(vector<Sequence> * seqs) {

	/* adds all sequences in the passed vector to the 'sequences'
	 * table in the db. */

	sqlite3 * db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		begin_transaction(db);

		// insert each keep sequence into the db, and set its db id in the keep_seqs vector
		for (int i = 0; i < seqs->size(); i++)
			seqs->at(i).set_sqlite_id(insert_sequence(db, seqs->at(i)));

		commit_transaction(db);
	} else
		cerr << "Could not add the sequence because the database " << name << " could not be opened" << endl;

	sqlite3_close(db);
}

int GeneDB::add_profile_alignment(string & profile_name, int child1_id, int child2_id, vector<Sequence> & sequences) {

	int new_profile_alignment_id = -1;

	sqlite3 *db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {

		begin_transaction(db);

		// create a new profile alignment record
		new_profile_alignment_id = insert_profile_alignment(db, profile_name, child1_id, child2_id);
		cout << "adding profile alignment: " << new_profile_alignment_id << endl;

		// insert the seqs into profile_alignment sequences
		vector<Sequence>::iterator s;
		for (s = sequences.begin(); s < sequences.end(); ++s) {

			int sequence_id = s->get_sqlite_id();
			string sequence = s->get_aligned_seq();
			insert_sequence_profile_map(db, sequence_id, new_profile_alignment_id, sequence);
		}
		commit_transaction(db);

	} else
		cerr << "Could not add the profile alignment because the database " << name << " could not be opened" << endl;

	sqlite3_close(db);

	return new_profile_alignment_id;

}

int GeneDB::add_profile_alignment(string & profile_name, vector<Sequence> & sequences) {

	/* just a wrapper for the full version of add_profile_alignment, for
	 * the case where there are no children (i.e. we are copying over
	 * an original alignment.
	 */

	int empty_child_id = 0;
	return add_profile_alignment(profile_name, empty_child_id, empty_child_id, sequences);
}

int GeneDB::add_empty_intermediate_profile(int child1_id, int child2_id) {

	/* this one is for compatibility with SQLiteProfiler, it calls this
	 * function and then loads the new alignment into the db itself. It might
	 * be better to deprecate this behavior and just have SQLiteProfiler call
	 * GeneDB to store the new profile alignments. this function is really
	 * just a wrapper for insert_profile_alignment()
	 */

	int new_profile_id = -1;
	string profile_name = to_string(0);

	sqlite3 * db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		begin_transaction(db);
		insert_profile_alignment(db, profile_name, child1_id, child2_id);
		commit_transaction(db);
		new_profile_id = sqlite3_last_insert_rowid(db);
	} else
		cerr << "Could not add new profile alignment because the database " << name << " could not be opened";

	sqlite3_close(db);

	return new_profile_id;
}

/*
 *
 * deprecated: used to be used by SQLiteProfiler, now it uses add_empty_intermediate_profile
 *
int GeneDB::add_profile_alignment(int child_id1, int child_id2) {

	// just creates a new record in the profile alignment table, and stores the ids
	// of its constituent child alignments

	// open db connection
	sqlite3 *conn2;
	int rc = sqlite3_open(name.c_str(), &conn2);
	char *zErrMsg = 0;
	sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
	if (rc != SQLITE_OK) {
		cerr << "Could not open database connection: " << *zErrMsg << endl;
		sqlite3_free(zErrMsg);
	}

	// build sql string
	string sql = "insert into profile_alignments (child1, child2, name) values (";
	sql += to_string(child_id1) + ",";
	sql += to_string(child_id2) + ",'";
	sql += to_string(0) + "');"; // not an ncbi (sas) - huh? (ceh)

	cout << sql << endl;

	// insert record
	rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
	int pid = sqlite3_last_insert_rowid(conn2);
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, &zErrMsg);
	if (rc != SQLITE_OK) {
		cerr << "Record insert failed: " << *zErrMsg << endl;
		sqlite3_free(zErrMsg);
	}

	// close the db connection
	sqlite3_close(conn2);

	// return the record's db id
	return pid;
}
*/

int GeneDB::add_alignment(string & filename, vector<Sequence> * dbseqs, vector<Sequence> * userseqs) {

	/* creates a new alignment record in the database 'alignments' table,
	 * and adds records to the sequence_alignment_map table to associate
	 * all the sequences in the dbseqs vector with the new alignment,
	 * then returns the new alignment's database id.
	 *
	 * if the new alignment record cannot be created, this should return
	 * -1 for the new alignment id, although it is likely that the program
	 * will abort due to an illegal operation before the return statement
	 * is reached.
	 */

	// TODO: add the reading of the actual alignment file
	int new_alignment_id = -1;

	sqlite3 * db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		cout << "adding alignment to database: " << name << endl;
		begin_transaction(db);

		// first create a new alignment record
		new_alignment_id = insert_alignment(db, filename, 0);

		// now add the dbseqs
		for (int i = 0; i < dbseqs->size(); i++) {
			int sequence_id = dbseqs->at(i).get_sqlite_id();
			string sequence = dbseqs->at(i).get_aligned_seq();
			insert_sequence_alignment_map(db, sequence_id, new_alignment_id, sequence);
		}

		// now the user seqs
		for (int i = 0; i < userseqs->size(); i++) {
			int sequence_id = userseqs->at(i).get_sqlite_id();
			string sequence = userseqs->at(i).get_aligned_seq();
			insert_sequence_alignment_map(db, sequence_id, new_alignment_id, sequence);
		}

		// clean up
		commit_transaction(db);
	} else
		cerr << "Could not create the alignment because the database connection could not be opened" << endl;

	sqlite3_close(db);

	return new_alignment_id;
}

void GeneDB::toggle_alignment_update(int alignment_id) {

	/* record that the alignment with the passed id has been updated */

	sqlite3 *db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		string sql = "update alignments set updated = 1 where id = ";
		sql += to_string(alignment_id) + ";";
		execute_simple_transaction(db, sql);
	} else
		cerr << "Could not record the alignment as updated because the database " << name << " could not be opened" << endl;

	sqlite3_close(db);
}

void GeneDB::remove_alignment(int alignment_id) {

	/* remove the alignment from the database */

	sqlite3 *db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		begin_transaction(db);

		string deleted_alignment_id = to_string(alignment_id);
		string sql = "delete from sequence_alignment_map where alignment_id = ";
		sql += deleted_alignment_id + ";";
		execute_query(db, sql);

		sql = "delete from alignments where id = ";
		sql += deleted_alignment_id + ";";
		execute_query(db, sql);

		commit_transaction(db);
	} else
		cerr << "Could not delete alignment " << alignment_id << " because the database " << name << " could not be opened" << endl;

	sqlite3_close(db);
}

void GeneDB::remove_profile_alignments() {

	sqlite3 *db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		begin_transaction(db);

		execute_query(db, "delete from sequence_profile_map;");
		execute_query(db, "delete from profile_alignments;");

		commit_transaction(db);
	} else
		cerr << "Could not clear the profile alignments because the database " << name << " could not be opened" << endl;

	sqlite3_close(db);
}

void GeneDB::get_first_profile_alignments(vector<string> & names) {

	/*
	 * these would be the childless alignments, the first ones
	 */

	// TODO: convert this to consistent style database interaction

	string sql = "select name from profile_alignments where child1 = 0 and child2 = 0;";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		names.push_back(to_string(query.getstr()));
	}
	query.free_result();
}

int GeneDB::get_alignment_id_by_name(string & alignname) {

	// TODO: convert this to consistent style database interaction

	string sql = "select id from alignments where alignname = '";
	sql += alignname + "';";

	Database conn(name);
	Query query(conn);
	int alignid;
	query.get_result(sql);
	while (query.fetch_row()) {
		alignid = atoi(to_string(query.getval()).c_str());
	}
	query.free_result();
	return alignid;
}

void GeneDB::remove_alignment_by_name(string alignname) {
	int alignid = get_alignment_id_by_name(alignname);
	remove_alignment(alignid);
}

void GeneDB::get_align_seqs(int alignid, vector<Sequence> & seqs) {

	/*
	 * seqs is returned with the id being the sqlite id and the seq being the
	 * current aligned seq
	 */

	// TODO: convert this to consistent style database interaction
	string alignids = to_string(alignid);
	string sql = "select sequence_id,sequence from sequence_alignment_map where alignment_id = ";
	sql += alignids + ";";

	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		string id;
		string seq;
		id = to_string(query.getval());
		seq = to_string(query.getstr());
		Sequence tseq(id, seq);
		seqs.push_back(tseq);
	}
	query.free_result();
}

void GeneDB::get_profile_align_seqs(int alignid, vector<Sequence> & seqs) {

	// TODO: convert this to consistent style database interaction

	string alignids = to_string(alignid);
	string sql = "select sequence_id,sequence from sequence_profile_map where profile_id = ";
	sql += alignids + ";";

	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		string id;
		string seq;
		id = to_string(query.getval());
		seq = to_string(query.getstr());
		Sequence tseq(id, seq);
		tseq.set_sqlite_id(atoi(id.c_str()));
		tseq.set_aligned_seq(seq);
		seqs.push_back(tseq);
	}
	query.free_result();
}

void GeneDB::get_align_seqs_unaligned(int alignid, vector<Sequence> & seqs) {

	/*
	 * seqs is returned with the id being the sqlite id and the seq being the
	 * original unaligned seq
	 */

	// TODO: convert this to consistent style database interaction
	string alignids = to_string(alignid);
	string sql =
			"select sequence_alignment_map.sequence_id,sequences.sequence from sequence_alignment_map,sequences where sequence_alignment_map.sequence_id=sequences.id and sequence_alignment_map.alignment_id = ";
	sql += alignids + ";";

	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		string id;
		string seq;
		id = to_string(query.getval());
		seq = to_string(query.getstr());
		Sequence tseq(id, seq);
		seqs.push_back(tseq);
	}
	query.free_result();
}

void GeneDB::get_align_seq_unaligned_fully_initialized(string alignname, vector<Sequence> & seqs) {

	/*
	 * seqs is returned with the id being the ncbi_id and the seq being the
	 * original unaligned seq and the other elements correctly defined
	 * the alignname is the name in the filename column
	 */

	// TODO: convert this to consistent style database interaction
	int alignid = get_alignment_id_by_name(alignname);
	string sql =
			"select sequence_alignment_map.sequence_id,sequences.ncbi_id,sequences.accession,sequences.name,sequences.sequence from sequence_alignment_map,sequences where sequence_alignment_map.sequence_id=sequences.id and sequence_alignment_map.alignment_id = ";
	sql += to_string(alignid) + ";";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		string id;
		string ncbi_id;
		string acc;
		string name;
		string seq;
		id = to_string(query.getval());
		ncbi_id = to_string(query.getval());
		acc = to_string(query.getstr());
		name = to_string(query.getstr());
		seq = to_string(query.getstr());
		Sequence tseq(ncbi_id, seq);
		if (name.find("user_") != string::npos) { //user id is username if user seq
			tseq.set_id(name);
		}
		tseq.set_ncbi_tax_id(ncbi_id);
		tseq.set_ncbi_gi_id(acc);
		tseq.set_name(name);
		tseq.set_sqlite_id(atoi(id.c_str()));
		seqs.push_back(tseq);
	}
	query.free_result();
}

void GeneDB::update_align_seqs(int alignment_id, vector<Sequence> & seqs) {

	/*
	 * uses the get_sqlite_id and the get_aligned_seq
	 */

	sqlite3 * db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		begin_transaction(db);

		for (int i = 0; i < seqs.size(); i++) {
			string sql = "update sequence_alignment_map set sequence = '";
			sql += seqs[i].get_aligned_seq() + "'";
			sql += " where alignment_id = ";
			sql += to_string(alignment_id) + " and sequence_id = ";
			sql += to_string(seqs[i].get_sqlite_id()) + ";";
			execute_query(db, sql);
		}
		commit_transaction(db);

	} else
		cerr << "Could not update the sequences because the database " << name << " could not be opened" << endl;

	sqlite3_close(db);
}

void GeneDB::update_profile_align_seqs(int alignid, vector<Sequence> & seqs) {

	/*
	 * uses the get_sqlite_id and the get_aligned_seq
	 */

	sqlite3 * db;
	if (sqlite3_open(name.c_str(), &db) == SQLITE_OK) {
		begin_transaction(db);

		for (int i = 0; i < seqs.size(); i++) {
			string sql = "update sequence_profile_map set sequence = '";
			sql += seqs[i].get_aligned_seq() + "'";
			sql += " where profile_id = ";
			sql += to_string(alignid) + " and sequence_id = ";
			sql += to_string(seqs[i].get_sqlite_id()) + ";";
			execute_query(db, sql);
			//cout << rc << " " << sql << endl;
		}

		commit_transaction(db);
	} else
		cerr << "Could not update the profile alignments because the database could not be opened" << endl;

	sqlite3_close(db);
}

void GeneDB::load_all_sequences_into(vector<Sequence> & seqs) {

	/*
	 * ncbi id will be the id and the comment
	 * name will be the name
	 * seq will be the original seq
	 * sqlite_id will be the sqlite_id
	 */

	// TODO: update to use consistent style database interaction

	// TODO: fix this so it is clear which fields are going where... it does not looking

	string sql = "select * from sequences;";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		string id;
		string ncbi_id;
		string acc;
		string name;
		string seq;
		id = to_string(query.getval());
		ncbi_id = to_string(query.getval());
		acc = to_string(query.getstr());
		name = to_string(query.getstr());
		seq = to_string(query.getstr());
		Sequence tseq(ncbi_id, seq);
		if (name.find("user_") != string::npos) { //user id is username if user seq
			tseq.set_id(name);
		}
		tseq.set_ncbi_tax_id(ncbi_id);
		tseq.set_ncbi_gi_id(acc);
		tseq.set_sqlite_id(atoi(id.c_str()));
		tseq.set_name(name);
		seqs.push_back(tseq);
	}
	query.free_result();
}

void GeneDB::load_orig_alignment_names_into(vector<string> & names) {

	// TODO: update to use consistent style database interaction

	string sql = "select alignname from alignments;";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		string name;
		name = to_string(query.getstr());
		names.push_back(name);
	}
	query.free_result();
}

void GeneDB::load_orig_alignment_nums_into(vector<int> & nums) {

	// TODO: update to use consistent style database interaction

	string sql = "select id from alignments;";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		int num;
		num = query.getval();
		nums.push_back(num);
	}
	query.free_result();
}

void GeneDB::load_first_profile_ids_into(vector<int> & nums) {

	/*
	 * typically used for updating and should be called after moving over updated
	 */

	// TODO: update to use consistent style database interaction

	string sql = "select id from profile_alignments where child1 = 0 and child2 = 0;";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while (query.fetch_row()) {
		int num;
		num = query.getval();
		nums.push_back(num);
	}
	query.free_result();
}

void GeneDB::load_sparse_sequences_from_original_alignment_into(vector<Sequence> & seqs, int alignment_id) {

	// open db connection
	Database conn(name);
	Query query(conn);

	// get the sequences for this alignment
	string sql = "select sequence_id, sequence from sequence_alignment_map where alignment_id = ";
	sql += to_string(alignment_id) + ";";
	query.get_result(sql);

	// for each row, populate a new sequence object and load it into the vector
	while (query.fetch_row()) {
		int dbid = query.getval();
		string seq = to_string(query.getstr());
		Sequence tseq(dbid,seq);
		seqs.push_back(tseq);
	}

	query.free_result();
}

void GeneDB::load_original_alignment_info_into(map<int, string> & the_map) {

	// open db connection
	Database conn(name);
	Query query(conn);

	// do query
	string sql = "select id, alignname from alignments;";
	query.get_result(sql);

	// load the results into the map
	while (query.fetch_row()) {
		int tid;
		string name;
		tid = query.getval();
		name = to_string(query.getstr());
		the_map[tid] = name;
	}

	query.free_result();
}

void GeneDB::migrate_original_alignments_and_load_new_profile_info_into(map<int, string> & profile_id_name_map) {

	// get info on all the original alignments
	map<int, string> original_alignments;
	load_original_alignment_info_into(original_alignments);

	map<int, string>::iterator a;
	for(a = original_alignments.begin(); a != original_alignments.end(); ++a) {

		int this_alignment_id = a->first;
		string this_alignment_name = a->second;

		// get the sequences associated this alignment
		vector<Sequence> sequences;
		load_sparse_sequences_from_original_alignment_into(sequences, this_alignment_id);

		// add a new profile alignment that is an exact duplicate of this original
		int new_profile_id = add_profile_alignment(this_alignment_name, sequences);

		// record the new profile info
		profile_id_name_map[this_alignment_id] = this_alignment_name;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      stopping here for today.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void GeneDB::copy_alignments_to_first_profiles_updated(map<int, string> & profile_id_name_map, vector<int>& updatedprofsnums) {
	string sql = "select id, alignname from alignments where updated = 1;";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	vector<int> aid;
	vector<string> names;
	while (query.fetch_row()) {
		int iid;
		string name;
		iid = query.getval();
		name = to_string(query.getstr());
		aid.push_back(iid);
		names.push_back(name);
	}
	query.free_result();
	//add the new alignment

	for (int i = 0; i < names.size(); i++) {
		sqlite3 *conn2;
		int rc = sqlite3_open(name.c_str(), &conn2);
		char *zErrMsg = 0;
		sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
		string sql = "insert into profile_alignments (child1,child2,name) values (0,0,'";
		sql += names[i] + "');";
		rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
		int pid = sqlite3_last_insert_rowid(conn2);
		sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
		updatedprofsnums.push_back(pid);

		sql = "select sequence_id,sequence from sequence_alignment_map where alignment_id = ";
		sql += to_string(aid[i]) + ";";
		Query query2(conn);
		query2.get_result(sql);
		vector<int> sid;
		vector<string> seqs;
		while (query2.fetch_row()) {
			sid.push_back(query2.getval());
			seqs.push_back(to_string(query2.getstr()));
		}
		query2.free_result();

		//insert the seqs into profile_alignment sequences
		sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
		for (int j = 0; j < sid.size(); j++) {
			sql = "insert into sequence_profile_map (sequence_id,profile_id,sequence) values (";
			sql += to_string(sid[j]) + ",";
			sql += to_string(pid) + ",'";
			sql += seqs[j] + "');";
			rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
		}
		sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
		sqlite3_close(conn2);
	}
	sql = "select id,name from profile_alignments where name != 0;";
	Query query3(conn);
	query3.get_result(sql);
	while (query3.fetch_row()) {
		int iid;
		string name;
		iid = query3.getval();
		name = to_string(query3.getstr());
		profile_id_name_map[iid] = name;
	}
	query3.free_result();

}

int GeneDB::get_deepest_profile_for_alignment(int alignid) {

	/* get the deepest alignment in the profile alignment table that includes the alignid.
	 * the deepest alignment should always be the one with the largest db id, because these
	 * ids increment as we progress down the guide tree (i.e. taxonomy). */

	// open db connection
	Database conn(name);
	Query query(conn);

	// retstack will hold alignment ids to be queried against
	stack<int> retstack;
	retstack.push(alignid);

	// this will hold the id of the deepest alignment, or NOT_YET_PROFILED if we find none
	int deepest = NOT_YET_PROFILED
	;

	while (!retstack.empty()) {
		int curid = retstack.top();
		retstack.pop();

		// get all profile alignments containing the alignment id'd by retstack.top()
		string sql = "select id from profile_alignments where ";
		sql += "child1 == " + to_string(curid) + " or ";
		sql += "child2 == " + to_string(curid) + ";";
//		cout << sql << endl;
		query.get_result(sql);

		// for every matching profile alignment
		while (query.fetch_row()) {
			int tid;
			string name;
			tid = query.getval();
			cout << "\tfound " << tid << endl;

			// if this profile alignment is the deepest yet; save it
			if (tid > deepest)
				deepest = tid;

			// add this profile to the retstack so we can check for
			// deeper alignments that contain it
			retstack.push(tid);
		}
		query.free_result();
	}
	cout << "deepest profile for " << alignid << ": " << deepest << endl;
	return deepest;
}


void GeneDB::write_profile_alignment_to_file(int alignid, string filename) {
	string sql = "select sequence_id,sequence from sequence_profile_map where profile_id = ";
	sql += to_string(alignid) + ";";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	vector<Sequence> seqs;
	while (query.fetch_row()) {
		string name;
		string seq;
		name = to_string(query.getval());
		seq = to_string(query.getstr());
//	cout << sql << " " << name << " " << seq << endl;
		Sequence tseq(name, seq);
		seqs.push_back(tseq);
	}
	query.free_result();
	FastaUtil fu;
	fu.writeFileFromVector(filename, seqs);
}

void GeneDB::write_profile_alignment_with_names_to_file(int alignid, string filename, bool ncbi) {
	string sql =
			"select sequence_profile_map.sequence_id,sequences.ncbi_id,sequences.name,sequence_profile_map.sequence from sequence_profile_map,sequences where sequence_profile_map.sequence_id=sequences.id and sequence_profile_map.profile_id = ";
	sql += to_string(alignid) + ";";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	vector<Sequence> seqs;
	while (query.fetch_row()) {
		string sid;
		string ncbiid;
		string name;
		string seq;
		sid = to_string(query.getval());
		ncbiid = to_string(query.getval());
		name = to_string(query.getstr());
		seq = to_string(query.getstr());
		if (ncbi && ncbiid != "0")
			name = ncbiid;
		Sequence tseq(name, seq);
		seqs.push_back(tseq);
	}
	query.free_result();
	FastaUtil fu;
	fu.writeFileFromVector(filename, seqs);
}

/*
 * assumes that the seq id is the sqlite id and the seq is the aligned seq
 */
void GeneDB::add_sequences_for_profile_alignment(int profilealignid, vector<Sequence> & seqs) {

	// set variables for db connection
	sqlite3 *conn2;
	int rc;
	char *zErrMsg = 0;

	// attempt to open the connection
	rc = sqlite3_open(name.c_str(), &conn2);
	rc = sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
	if (rc != SQLITE_OK) {
		cerr << "SQL error: %s\n" << zErrMsg;
		sqlite3_free(zErrMsg);
	}

	for (int i = 0; i < seqs.size(); i++) { // insert the sequences
		string sql = "insert into sequence_profile_map (sequence_id,profile_id,sequence) values (";
		sql += to_string(seqs[i].get_id()) + ",";
		sql += to_string(profilealignid) + ",'";
		sql += to_string(seqs[i].get_sequence()) + "');"; //not an ncbi

		// execute the query
		rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, &zErrMsg);
		if (rc != SQLITE_OK) {
			cerr << "SQL error: %s\n" << zErrMsg;
			sqlite3_free(zErrMsg);
		}
	}
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, &zErrMsg);
	if (rc != SQLITE_OK) {
		cerr << "SQL error: %s\n" << zErrMsg;
		sqlite3_free(zErrMsg);
	}
	sqlite3_close(conn2);
}

void GeneDB::get_updated_profs_names_delete_old(vector<string> & updatedprofs, vector<int> & updatedprofsnums, vector<int>&notupdatedprofnums) {
	//get the updated names
	string sql =
			"select alignments.id,alignments.alignname,profile_alignments.id from alignments,profile_alignments where alignments.alignname=profile_alignments.name and alignments.updated = 1;";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	stack<int> retstack;
	vector<int> todelete;
	while (query.fetch_row()) {
		int sid;
		string name;
		int aid;
		sid = query.getval();
		name = to_string(query.getstr());
		aid = query.getval();
		updatedprofs.push_back(name);
		updatedprofsnums.push_back(aid);
		retstack.push(aid);
		todelete.push_back(int(aid));
	}
	query.free_result();
	sql = "select id,name from profile_alignments where child1 = 0 and child2 = 0;";
	Query queryt(conn);
	queryt.get_result(sql);
	while (queryt.fetch_row()) {
		int sid;
		string name;
		sid = queryt.getval();
		name = to_string(queryt.getstr());
//	cout << "notupdatedsearch: " << sid << " " << name << endl;
		if (count(updatedprofs.begin(), updatedprofs.end(), name) == 0) {
			notupdatedprofnums.push_back(int(sid));
//	    cout << "sid: " << sid << endl;
		}
	}
	queryt.free_result();

	//get all the profiles involving these names
	int finalret = -1;
	while (!retstack.empty()) {
		int curid = retstack.top();
		retstack.pop();
		string sql = "select id from profile_alignments where child1 = ";
		sql += to_string(curid) + " or child2 = ";
		sql += to_string(curid) + ";";
		Database conn(name);
		Query query(conn);
		query.get_result(sql);
		while (query.fetch_row()) {
			int tid;
			string name;
			tid = query.getval();
			if (count(todelete.begin(), todelete.end(), int(tid)) == 0)
				todelete.push_back(tid);
			retstack.push(tid);
		}
		query.free_result();
	}

	//remove all the profiles that involve the names
	//TODO: could get the process
	for (int i = 0; i < todelete.size(); i++) {
		sqlite3 *conn2;
		int rc = sqlite3_open(name.c_str(), &conn2);
		char *zErrMsg = 0;
		sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
		sql = "delete from sequence_profile_map where profile_id =";
		sql += to_string(todelete[i]) + ";";
		rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
		sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);

		sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
		sql = "delete from profile_alignments where id =";
		sql += to_string(todelete[i]) + ";";
		rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
		sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
		sqlite3_close(conn2);
	}
	load_first_profile_ids_into(updatedprofsnums);
}

void GeneDB::toggle_updated_all_off() {
	sqlite3 *conn;
	int rc = sqlite3_open(name.c_str(), &conn);
	char *zErrMsg = 0;

	sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
	string sql = "update alignments set updated = 0;";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
	sqlite3_close(conn);
}
