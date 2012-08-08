#!/usr/bin/env python

import os, sys, re, sqlite3, subprocess
from decimal import *
from copy import deepcopy

# using Decimal class for high precision
getcontext().prec = 100

# default value for some parameters
BIGNUMBER = 10000;

class Sequence():
    # just a container for convenient access
    def __init__(self, info):
        self.gi = info[0]
        self.taxon_name = info[1]
        self.description = info[2]
        self.sequence = info[3]

def get_any_depth_children_by_rank(db, parent, rank):
    # returns a list of tuples describing all children
    # of parent of specified rank 

    con = sqlite3.connect(db)
    cur = con.cursor()

    # get parent taxon
    sql = "SELECT ncbi_id, left_value, right_value FROM taxonomy WHERE name == ? AND name_class = 'scientific name';"
    cur.execute(sql, (parent, ))
    r = cur.fetchone()
   
    ncbi_id = r[0]
    parent_lval = r[1]
    parent_rval = r[2]
    print "Found one taxon matching '" + parent + "' with ncbi_id: " + str(ncbi_id)

    # get all child taxa of specified rank
    cur.execute("SELECT name, ncbi_id, left_value, right_value FROM taxonomy WHERE " \
                    "left_value > ? AND right_value < ? AND node_rank == ? AND " \
                    "name_class == 'scientific name';", (parent_lval, parent_rval, rank))
    records = cur.fetchall()

    # record children
    children = list()
    for r in records:
        children.append({"name": r[0], "ncbi_id": r[1], "lval": r[2], "rval": r[3]})

    con.close()
    return children

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "usage: find_guide_sequences.py <cfgfile>"
        sys.exit(0)
 
    # default parameter values
    pars = {"minexpect": 0, \
            "maxexpect": Decimal(BIGNUMBER), \
            "minscore": 0, \
            "maxscore": Decimal(BIGNUMBER), \
            "minlength": 200, \
            "maxlength": 5000, \
            "nametargets": False}

    # read parameters from config file
    cfgfile = open(sys.argv[1],"r")
    for line in cfgfile:
        
        # parse line
        toks = [t.strip() for t in line.split("=")]
        parname = toks[0]
        if len(toks) > 1:
            val = toks[1]

        # for raw text parameters
        if parname in ["db","clade","targetrank","seedfile"]:
            pars[parname] = val

        # for numeric parameters
        elif parname in ["maxexpect","minexpect","maxscore","minscore","minlength"]:
            pars[parname] = float(val)

        # needs special formatting
        elif parname == "search":
            pars[parname] = "%" + val.strip("\"").strip("'") + "%"

        # boolean parameters
        elif parname in ["nametargets",]:
            pars[parname] = True
    
    required = ["db", "clade", "targetrank", "seedfile", "search"]
    if not all(a in pars.keys() for a in required):
        print "config file is missing at least one required parameter (" + ", ".join(required) + ")"
        sys.exit(0)

    # set up a temp directory to store fasta files for bl2seq
    tempdir_name = "guidemaker_temp"    
    try:
        os.mkdir(tempdir_name)
    except OSError:
        pass
    tempdir = os.path.basename(tempdir_name) + "/"

    # read sequences from seedfile
    seedfile = open(pars["seedfile"],'r')
    seeds = dict()
    n_seeds = 0
    eof = False
    
    while True:
        line = seedfile.readline()

        if len(line) == 0:
            # end of file; attempt to write last seed before break
            writelast = True
            eof = True
        elif line[0] == '>':
            # hit a new sequence; finish the last one and move on
            writelast = True
            startnew = True
        else:
            # not a new sequence or eof, so append to current sequence
            writelast = False
            startnew = False

        if writelast and "seed_info" in locals().keys():
            # if there is a previous seed waiting to be written, do it 
            last_seed = Sequence(seed_info)
            this_seedfile = tempdir + "seed." + \
                            "_".join((pars["clade"],"by",pars["targetrank"],pars["search"].strip("%").strip())) + \
                            "." + str(n_seeds) + ".fasta"
            this_seed_file = open(this_seedfile, "w")
            this_seed_file.write(">" + last_seed.taxon_name + "\n" + last_seed.sequence)
            this_seed_file.close()
            seeds[n_seeds] = last_seed.taxon_name
            n_seeds += 1

        if startnew:
            # seed_info[0] and [2] are empty because we don't care
            # about gi numbers and descriptions for seeds
            seed_seq_name = line[1:].strip()
            seed_info = [0, seed_seq_name, "", ""]
        else:
            # not a new seq or the eof, so append
            seed_info[3] += line.strip()
            
        if eof:
            break

    # container to hold best guide seq sets: n_seeds empty lists, indexed by seed number
    # need to deepcopy rows so they don't all reference the same object
    empty_rows = map(deepcopy, [{"score": [], "expect": []},] * n_seeds)
    proposed_sets = dict(zip(range(n_seeds),empty_rows))

    # find target taxa
    targets = get_any_depth_children_by_rank(pars["db"], pars["clade"], pars["targetrank"])

    con = sqlite3.connect(pars["db"])
    cur = con.cursor()
    
    # get count of candidate seqs for this parent
    # should make this more elegant so we're not searching for parent twice
    cur.execute("SELECT left_value, right_value FROM taxonomy WHERE " \
                "name == ? AND name_class = 'scientific name';", (pars["clade"], ))
    parent_lval, parent_rval = cur.fetchone()
    cur.execute("SELECT count(*) FROM sequence WHERE ncbi_id IN " \
                "(SELECT ncbi_id FROM taxonomy WHERE left_value > ? AND " \
                "right_value < ?) AND description LIKE ?;", \
                 (parent_lval, parent_rval, pars["search"]))
    n_rows_all = cur.fetchone()[0]
    print "Will be evaluating " + str(n_rows_all) + " possible guide sequences from " + pars["clade"]

    i = 0
    for target in targets:

        # to record the best score and expect value for each seed
        bestscore_for_seed = dict(zip(range(n_seeds),[0,]*n_seeds))
        # deepcopy rows so they don't all reference same object
        bestexpect_for_seed = dict(zip(range(n_seeds),map(deepcopy,[Decimal(BIGNUMBER),]*n_seeds)))

        # to record the gi number of the candidate with the best score/expect for each seed
        bestscore_gi_for_seed = dict(zip(range(n_seeds),[None,]*n_seeds))
        bestexpect_gi_for_seed = dict(zip(range(n_seeds),[None,]*n_seeds))

        if pars["nametargets"]:
            # get count of candidate seqs for this target
            cur.execute("SELECT count(*) FROM sequence WHERE ncbi_id IN " \
                        "(SELECT ncbi_id FROM taxonomy WHERE left_value > ? AND " \
                        "right_value < ?) AND description LIKE ?;", \
                         (target["lval"], target["rval"], pars["search"],))
            n_rows_target = cur.fetchone()[0]
            if i > 0:
                print "\n"
            print "Found " + str(n_rows_target) + " matching rows for family " + target["name"]
            next_target_i = i + n_rows_target + 1
            next_target_message = ("(next target at " + str(next_target_i) + ")").rjust(30)
        else:
            next_target_message = ""

        # get all candidate seqs for this target
        cur.execute("SELECT identifier, name, description, seq FROM " \
                    "sequence INNER JOIN taxonomy ON sequence.ncbi_id == taxonomy.ncbi_id WHERE " \
                    "name_class == 'scientific name' AND left_value > ? AND right_value < ? AND " \
                    "description LIKE ?;", \
                     (target["lval"], target["rval"], pars["search"]))
        result = cur.fetchall()

        # for each candidate sequence
        for seqinfo in result:
            i += 1

            # use Sequence object to parse query result for this candidate
            candidate = Sequence(seqinfo)

            # skip seqs not within range of acceptable lengths
            if len(candidate.sequence) < pars["minlength"] or len(candidate.sequence) > pars["maxlength"]:
                continue

            # user feedback
            sys.stdout.write("\r{} / {} total {} ".format(i, n_rows_all, next_target_message))
            sys.stdout.flush()

            # write the candidate to a temp file
            tempfile = open(tempdir + candidate.gi + ".fasta", "w")
            tempfile.write(">" + candidate.description + "\n" + candidate.sequence)
            tempfile.close()
            
            # use blastn to compare to each seed
            for s in range(n_seeds):
                bl2seq_args = ["bl2seq", \
                    "-i", tempdir + "seed." +  \
                            "_".join((pars["clade"],"by",pars["targetrank"],pars["search"].strip("%").strip())) + \
                            "." + str(s) + ".fasta",
                    "-j", tempdir + candidate.gi + ".fasta", \
                    "-p", "blastn"]
                p1 = subprocess.Popen(bl2seq_args, stdout=subprocess.PIPE)
                bl2seq_out = p1.communicate()[0]

                # extract score info from the bl2seq output
                candidate.score = 0
                candidate.expect = Decimal(BIGNUMBER)

                # skip results returning no hits
                testlines = bl2seq_out.split("\n",4)
                if (testlines[3].strip() + " ")[0] == "*":
                    continue
                # this might be faster, i have not tested it
#                if not bl2seq_out.find("***** No hits found ******") >= 0:
#                    continue

                # find score info for all hits in the bl2seq
                score_sets = re.finditer(r"Score \= ([0-9.]+) bits \([0-9.]+\)\, Expect \= ([e0-9\-.]+)",bl2seq_out)
                for scores in score_sets:
                    c, e = [float(a) for a in scores.groups()]

                    # sum the score values for all hits of this candidate to the current seed
                    candidate.score += c

                    if Decimal(e) < candidate.expect:
                        # find the best expected value for this sequence
                        candidate.expect = Decimal(e)
                
                # record the score if it is the best yet
                if candidate.score > bestscore_for_seed[s]:
                    # unless it is still not better than the min allowed score
                    if candidate.score > pars["minscore"]:
                        bestscore_for_seed[s] = candidate.score
                        bestscore_gi_for_seed[s] = candidate.gi
                
                # record the expect value if it is the best yet
                if candidate.expect < bestexpect_for_seed[s]:
                    # unless it is still not better than the max allowed expect
                    if candidate.expect < pars["maxexpect"]:
                        bestexpect_for_seed[s] = candidate.expect
                        bestexpect_gi_for_seed[s] = candidate.gi

                # still need to add mechanism that stops searching once we have hit a sequence
                # with a score > maxscore or an expect < minexpect

        # store best gis in proposed_sets
        for seed in range(n_seeds):
            proposed_sets[seed]["score"].append(bestscore_gi_for_seed[seed])
            proposed_sets[seed]["expect"].append(bestexpect_gi_for_seed[seed])

    print ""

    outfile = open(pars["clade"] +"_by_" + pars["targetrank"] + ".guideseqs", "w")
    for s in range(n_seeds):
        outfile.write("seed " + seeds[s] + " by score:\n")
        for gi in proposed_sets[s]["score"]:
            if gi != None:
                outfile.write(gi+"\n")
        outfile.write("\n")

        outfile.write("seed " + seeds[s] + " by expect:\n")
        for gi in proposed_sets[s]["expect"]:
            if gi != None:
                outfile.write(gi+"\n")
        outfile.write("\n")
    
    outfile.close()
