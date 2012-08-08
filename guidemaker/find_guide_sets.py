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
        print "usage: find_guide_sequences.py <path to .guide file or dir with .guide files>"
        sys.exit(0)
 
    # default parameter values
    pars = {"minexpect": 0, \
            "maxexpect": Decimal(BIGNUMBER), \
            "minscore": 0, \
            "maxscore": Decimal(BIGNUMBER), \
            "minlength": 200, \
            "maxlength": 5000, \
            "nametargets": False}

    # read config file(s)
    cfgfiles = []
    if not os.path.isdir(sys.argv[1]):
        workingdir = ""
        cfgfiles.append(sys.argv[1])
    else:
        workingdir = sys.argv[1].strip("/") + "/"
        cfgfiles = [fname for fname in os.listdir(workingdir) if fname[-6:] == ".guide"]

    # run over each config file
    for cfgfname in cfgfiles:

        cfgfile = open(workingdir+cfgfname,"r")
        for line in cfgfile:
            
            # parse line
            toks = [t.strip() for t in line.split("=")]
            parname = toks[0]
            if len(toks) > 1:
                val = toks[1]

            # for raw text parameters
            if parname in ["db","targetrank","seedfile"]:
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
            
            # additional processing
            if parname == "clades":
                pars[parname] = [v.strip() for v in val.split(",")]
        
        required = ["db", "clades", "targetrank", "seedfile", "search"]
        if not all(a in pars.keys() for a in required):
            print "config file is missing at least one required parameter (" + ", ".join(required) + ")"
            sys.exit(0)

        idtag = "_".join(("_".join(pars["clades"]),"by",pars["targetrank"],pars["search"].strip("%").strip()))
        print "output will be labeled as " + idtag

        # set up a temp directory to store fasta files for bl2seq
        tempdir_name = "guidemaker_temp"    
        try:
            os.mkdir(tempdir_name)
        except OSError:
            pass
        tempdir = os.path.basename(tempdir_name) + "/"

        # open seedfile
        seedfile = open(pars["seedfile"],'r')
        seeds = dict()
        n_seeds = 0
        eof = False
            
        # read seed sequences
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
                this_seedfile = tempdir + "seed." + str(n_seeds) + "." + idtag + ".fasta"
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
        proposed_set_for_seed = dict(zip(range(n_seeds),map(deepcopy,[[],]*n_seeds)))

        # find target taxa
        targets = []
        for clade in pars["clades"]:
            targets += get_any_depth_children_by_rank(pars["db"], clade, pars["targetrank"])

        con = sqlite3.connect(pars["db"])
        cur = con.cursor()
        
        # get count of candidate seqs for this parent
        # should make this more elegant so we're not searching for parent twice
        n_rows_all = 0
        for clade in pars["clades"]:
            cur.execute("SELECT left_value, right_value FROM taxonomy WHERE " \
                        "name == ? AND name_class = 'scientific name';", (clade, ))
            parent_lval, parent_rval = cur.fetchone()
            cur.execute("SELECT count(*) FROM sequence WHERE ncbi_id IN " \
                        "(SELECT ncbi_id FROM taxonomy WHERE left_value > ? AND " \
                        "right_value < ?) AND description LIKE ?;", \
                         (parent_lval, parent_rval, pars["search"]))
            n_rows_all += cur.fetchone()[0]
        print "Will be evaluating " + str(n_rows_all) + " possible guide sequences from " + ", ".join(pars["clades"])

        i = 0
        for target in targets:

            # to record the best expect value for each seed
            # deepcopy rows so they don't all reference same object
            bestexpect_for_seed = dict(zip(range(n_seeds),map(deepcopy,[Decimal(BIGNUMBER),]*n_seeds)))

            # to record the gi number of the candidate with the best expect for each seed
            bestexpect_gi_for_seed = dict(zip(range(n_seeds),[None,]*n_seeds))

            if pars["nametargets"]:
                # report count of candidate seqs for this target
                cur.execute("SELECT count(*) FROM sequence WHERE ncbi_id IN " \
                            "(SELECT ncbi_id FROM taxonomy WHERE left_value > ? AND " \
                            "right_value < ?) AND description LIKE ?;", \
                             (target["lval"], target["rval"], pars["search"],))
                n_rows_target = cur.fetchone()[0]
                if i > 0:
                    print "\n"
                print "Found " + str(n_rows_target) + " matching rows for target " + target["name"]
                next_target_i = i + n_rows_target + 1
                next_target_message = ("(next target at " + str(next_target_i) + ")").rjust(30)
            else:
                next_target_message = ""

            # retrieve candidate seqs for this target
            cur.execute("SELECT identifier, name, description, seq FROM " \
                        "sequence INNER JOIN taxonomy ON sequence.ncbi_id == taxonomy.ncbi_id WHERE " \
                        "name_class == 'scientific name' AND left_value > ? AND right_value < ? AND " \
                        "description LIKE ?;", \
                         (target["lval"], target["rval"], pars["search"]))
            result = cur.fetchall()

            # for each candidate sequence
            for seqinfo in result:

                # user feedback
                i += 1
                sys.stdout.write("\r{} / {} total {} ".format(i, n_rows_all, next_target_message))
                sys.stdout.flush()

                # use Sequence object to parse query result for this candidate
                candidate = Sequence(seqinfo)

                # skip seqs not within range of acceptable lengths
                if len(candidate.sequence) < pars["minlength"] or len(candidate.sequence) > pars["maxlength"]:
                    continue

                # write the candidate to a temp file
                tempfile = open(tempdir + candidate.gi + ".fasta", "w")
                tempfile.write(">" + candidate.description + "\n" + candidate.sequence)
                tempfile.close()
                
                # use blastn to compare to each seed
                for s in range(n_seeds):
                    bl2seq_args = ["bl2seq", \
                        "-i", tempdir + "seed." + str(s) + "." + idtag + ".fasta",
                        "-j", tempdir + candidate.gi + ".fasta", \
                        "-p", "blastn"]
                    p1 = subprocess.Popen(bl2seq_args, stdout=subprocess.PIPE)
                    bl2seq_out = p1.communicate()[0]

                    # extract score info from the bl2seq output
                    candidate.expect = Decimal(BIGNUMBER)

                    # skip results returning no hits
                    testlines = bl2seq_out.split("\n",4)
                    if (testlines[3].strip() + " ")[0] == "*":
                        continue
                    # this might be faster, i have not tested it
    #                if not bl2seq_out.find("***** No hits found ******") >= 0:
    #                    continue

                    # find the best expect value for all hits to this sequence
                    expect_vals = re.finditer(r"Score \= [0-9.]+ bits \([0-9.]+\)\, Expect \= ([e0-9\-.]+)",bl2seq_out)
                    for e in [float(match.groups()[0]) for match in expect_vals]:
                        if Decimal(e) < candidate.expect:
                            candidate.expect = Decimal(e)
                        
                    # record the expect value if it is the best yet
                    if candidate.expect < bestexpect_for_seed[s]:
                        # unless it is still not better than the max allowed expect
                        if candidate.expect < pars["maxexpect"]:
                            bestexpect_for_seed[s] = candidate.expect
                            bestexpect_gi_for_seed[s] = candidate.gi

                    ### still need to add mechanism that stops searching for matches
                    ### to a particular seed once we have hit a sequence for that seed
                    ### with an expect < minexpect. we could just maintain a variable
                    ### somewhere, indexed by seed, that would hold a true value if a
                    ### best match had been found and a false if not

            # add best candidate gi for this target to proposed sets
            for seed in range(n_seeds):
                proposed_set_for_seed[seed].append(bestexpect_gi_for_seed[seed])

        # makes ouput look nicer
        print ""

        ### should also have option to identify highest-scoring set of guide
        ### seqs by x-blasting all seqs within each set

        try:
            os.mkdir("keepfiles")
        except OSError:
            pass

        # write proposed sets to file
        guidesetfile = open("keepfiles/" + idtag + ".guideseqs", "w")
        for s in range(n_seeds):
            guidesetfile.write("seed " + seeds[s].split("|",2)[1] + " by expect:\n")
            for gi in proposed_set_for_seed[s]:
                if gi != None:
                    guidesetfile.write(gi+"\n")
            guidesetfile.write("\n")
        guidesetfile.close()
        
        for s in range(n_seeds):
            keepfile = open("keepfiles/" + idtag + "." + seeds[s].split("|",2)[1] + ".keep", "w")
            for gi in proposed_set_for_seed[s]:
                if gi != None:
                    cur.execute("SELECT description, seq from sequence where identifier == ?;", (gi,))
                    desc, seq = cur.fetchone()
                    keepfile.write(">gi|"+gi+"||| "+desc+"\n"+seq+"\n")
            keepfile.close()
        con.close()
