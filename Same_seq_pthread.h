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
 *  Same_seq_pthread.h
 */

#ifndef _SAME_SEQ_PTHREAD_H_
#define _SAME_SEQ_PTHREAD_H_

#include <string>
#include <vector>

using namespace std;

#include <Seq/Sequence.h>
#include <Seq/containers>
#include <Seq/ioseq>
#include <Seq/alphabets>
using namespace bpp;

#include "DBSeq.h"
#include "utils.h"

struct thread_data{
	int thread_id;
	vector<DBSeq> seqs;
	vector<DBSeq> keep_seqs;
	vector<bool> keep_rc;
	int reports;
	double coverage;
	double identity;
	OrderedSequenceContainer * known_seqs;
};

vector<double> get_blast_score_and_rc_cstyle(Sequence inseq1, DBSeq inseq2, bool * prevcomp, int thread);

void * Same_seq_pthread_go(void *threadarg);

#endif
