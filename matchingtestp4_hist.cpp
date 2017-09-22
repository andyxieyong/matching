/*
 * matchingtestp4_hist.cpp
 * 
 * match transportation and CDR data
 * use the crecords class for storing CDR data and searching
 * read LTA records from text file, iterate over it, store possible / impossible matches in hash sets
 * 	(will probably use a lot of memory)
 * 
 * parallel implementation
 * 
 * improved version, use two crecords classes (one for CDR and one for LTA), perform matching in both directions (if needed)
 * 
 * further improved version, use copies of both dataset sorted by uids, first generate possible matches then check for impossible matches among these
 *  + extra option to test if old and new versions of functions give the same result
 * 
 * read data from binary storage files
 * 
 * new possibilities:
 * 	keep track of all points, allow only one match per point in both datasets
 * 
 * create statistics of only possible or with impossible matches
 * 
 * output only aggregate histograms possibly grouped by the number of records of the users
 * 	(similarly to the random match program)
 * 
 * Copyright 2016 Kondor Dániel <dkondor@mit.edu>
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following disclaimer
 *   in the documentation and/or other materials provided with the
 *   distribution.
 * * Neither the name of the  nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * 
 */

#ifndef USECEREAL
#define USECEREAL
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <unordered_map> // use C++11 STL hashmaps instead of google::dense_hash_*
#include <unordered_set>
#include <set>
#include <signal.h>
#include <unistd.h>
#include <algorithm> // std::sort for determining the top matches
#include <pthread.h>
#include "cdrrecords.h"
#include "popen_noshell.h"



struct sbounds {
		uint32_t ucmin;
		uint32_t ucmax;
		uint32_t utmin;
		uint32_t utmax;
		sbounds():ucmin(0),ucmax(0),utmin(0),utmax(0) {  }
		sbounds(uint32_t _ucmin, uint32_t _ucmax, uint32_t _utmin, uint32_t _utmax):ucmin(_ucmin),ucmax(_ucmax),utmin(_utmin),utmax(_utmax) {  }
		
		void do_copy_basic(const sbounds& b) {
			ucmin = b.ucmin;
			ucmax = b.ucmax;
			utmin = b.utmin;
			utmax = b.utmax;
		}
		
		sbounds(const sbounds& b) {
			do_copy_basic(b);
		}
		sbounds(const sbounds&& b) {
			do_copy_basic(b);
		}
		sbounds& operator=(const sbounds& b) {
			do_copy_basic(b);
			return *this;
		}
		
		// read one unsigned int number from the supplied FILE* stream, skipping spaces and tabs but not newlines
		static inline bool readnum(FILE* fin, uint32_t* r) {
			int a;
			do { a = getc(fin); } while(a == ' ' || a == '\t');
			if(isdigit(a)) {
				ungetc(a,fin);
				return ( fscanf(fin,"%u",r) == 1 );
			}
			else return false;
		}
		inline bool readbounds(FILE* fin) {
			if( readnum(fin,&ucmin) ) if( readnum(fin,&ucmax) ) if( readnum(fin,&utmin) ) if( readnum(fin,&utmax) ) return true;
			return false;
		}
		inline void printbounds(FILE* fout) {
			fprintf(fout,"%u\t%u\t%u\t%u",ucmin,ucmax,utmin,utmax);
		}
		
		// compare: first ucmin, then utmin
		// important: it is assumed that bounds do not overlap!
		// an exception is thrown if this is not the case
		bool operator < (const sbounds& s) const {
			if(ucmin < s.ucmin) {
				if(ucmax >= s.ucmin) throw new std::runtime_error("sbounds::operator <(): overlapping bounds supplied!\n");
				else return true;
			}
			// ucmin >= s.ucmin
			if(ucmin == s.ucmin) {
				if(ucmax != s.ucmax) throw new std::runtime_error("sbounds::operator <(): overlapping bounds supplied!\n");
				if(utmin < s.utmin) {
					if(utmax >= s.utmin) throw new std::runtime_error("sbounds::operator <(): overlapping bounds supplied!\n");
					else return true;
				}
				if(utmin == s.utmin) throw new std::runtime_error("sbounds::operator <(): non-unique bounds!\n");
				// utmin > s.utmin
				if(utmin <= s.utmax) throw new std::runtime_error("sbounds::operator <(): overlapping bounds supplied!\n");
				return false;
			}
			// ucmin > s.ucmin
			if(ucmin <= s.ucmax) throw new std::runtime_error("sbounds::operator <(): overlapping bounds supplied!\n");
			return false;
		}
		
		// compare with actual record numbers of users (to choose an appropriate bin)
		// see if the supplied user numbers actually fall in this bin
		bool contains(uint32_t nc, uint32_t nt) const {
			if(nc >= ucmin && nc <= ucmax && nt >= utmin && nt <= utmax) return true;
			else return false;
		}
	};
	
	// compare with actual record numbers of users (to choose an appropriate bin)
	// use these with std::lower_bound
	bool operator < (const sbounds& s, const std::pair<uint32_t,uint32_t>& n) {
		if(s.ucmax < n.first) return true;
		if(s.ucmin <= n.first) if(s.utmax < n.second) return true;
		return false;
	}
	
	bool operator < (const std::pair<uint32_t,uint32_t>& n, const sbounds& s) {
		if(s.ucmin > n.first) return true;
		if(s.ucmax >= n.first) if(s.utmin > n.second) return true;
		return false;
	}
	
	
	

// shared data by the threads (records to process) -- note: this could be implemented as a derived class of crecords
struct rdata {
	//distinct user ids and pointers in the above vectors to the start of the data associated with them
	std::vector<int64_t> uids_uniq;
	
	//number of users processed (a thread should choose it's next user to process as uids_uniq[nusers]
	unsigned int nusers;
	//~ pthread_mutex_t nusers_mutex; //mutex to use to choose the next user to process
	
	//~ bool out_impossible_brief; //only write one record per user to out_impossible, the number of matches found
	//~ bool quick; //do not keep track of the number of matches for impossible combinations
	unsigned int dt1; // time window to use for walking
	unsigned int dt2; // time window to use for transit
	
	// 1. store data sorted by cell_id and timestamp for searching for possible matches
	crecords* cr; //store cdr records here (only read by the threads, not modified)
	// 2. store data to search for sorted by user_id and timestamp + use the same data for checking impossible matches
	crecords* tr; //store transportation data here (readonly too)
	// 3. store the CDR data again sorted by uid and timestamp to check for impossible matches
	crecords* cru;
	bool searchreverse; // if true, cr and cru store transportation records and tr store CDR records; the compare_users*() functions need to be called on cru
	
	unsigned int out_p_max; //output the first this many possible matches only
	
	bool skiperrors; //if true, errors from search_records_possible() (i.e. missing place ID) do not result in an aborted run
		//when searching in reverse (searching for CDR records among the LTA data), there can be cells with no stops in or around them,
		//that is not an error
	bool statall; //if true, all temporal matches are counted (merged the functionality from matchingtestp4_all.cpp); cr needs to be sorted by time
		// in this case, temporal matches are output as well
	bool onematch; //if true, each point can result in only one match also in the searched dataset
	bool error; // signal an error during processing
	
	uint32_t tmin; // if > 0, limit searches to this range
	uint32_t* tmax; // if ntmax > 0 use several search ranges
	int ntmax;
	
	// limit the number of running threads and / or inform the main thread about the number of running threads
	// (need to be set to the number of threads before start, thread decrement these when exiting)
	unsigned int nthreads;
	unsigned int running_threads;
	unsigned int tid;
	unsigned int out_progress_uid; // output progress after every this many uids
	
	// vector of bounds to consider in the result histograms (readonly)
	std::vector<sbounds>* vbounds;
	
	// result (output; only created in the end by summing up output of every thread)
	pthread_mutex_t res_mutex; // mutex for writing to the result histograms
	std::vector<int64_t>* res_possible_matches; // vectors with possible matches:
		// res_possible_matches[ vbounds.size() * i + j ][ k ] contains the number of pairs with k+1 matches in time period [tmin,tmax[i]) and
		// between users in vbounds[j] part
		// in total it has to be allocated to contain ntmax*vbounds.size() empty vectors in the beginning
	std::vector<int64_t>* res_temporal_matches; // vectors with temporal matches, laid out similarly to the previous one (only used if statall == true)
};


// run this function on multiple threads, do the actual processing here
void* process_records(void* p0) {
	
	// set thread signal mask to ignore SIGINT (that should be handled in the main thread)
	sigset_t set;
	sigemptyset(&set);
	sigaddset(&set, SIGQUIT);
	sigaddset(&set, SIGINT);
	pthread_sigmask(SIG_BLOCK, &set, 0);
	
	rdata* rd = (rdata*)p0;
	
	// sets used during normal processing
	std::unordered_set<int64_t> s_candidates; //use this to store all / any matches which are later checked -- if all temporal matches are counted
	
	const size_t nres = (rd->ntmax) * (rd->vbounds->size());
	size_t nemax = 0; // size of the vectors below
	std::vector<uint64_t>* h_possible_matches   = new std::vector<uint64_t>[nres]; //used as temporary storage to create histograms of the number of matches
	std::vector<uint64_t>* h_impossible_matches = new std::vector<uint64_t>[nres]; //note: these are resized dynamically as needed
	unsigned int* matches_possible = 0;
	unsigned int* matches_temporal = 0;
	if(rd->ntmax > 1) {
		matches_possible = new unsigned int[rd->ntmax];
		matches_temporal = new unsigned int[rd->ntmax];
	}
	
	// assign thread ID
	unsigned int tid1 = __atomic_fetch_add(&(rd->tid),1,__ATOMIC_ACQ_REL);
	bool running = true;
	
	const std::vector<sbounds>& vbounds = *(rd->vbounds); // easier management
	
	while(1) {
		// check for error or quit request -- optionally stop processing if requested
		if( rd->error || rd->nthreads == 0 ) break;
		if( tid1 >= rd->nthreads ) {
			running = false;
			unsigned int nt2 = __atomic_fetch_sub(&(rd->running_threads),1,__ATOMIC_ACQ_REL);
			if(nt2 == 0) { // subtracting 1 from unsigned 0 is undefined behavior, this is an error then
				rd->running_threads = 0;
				rd->error = true;
				break;
			}
			
			do sleep(5);
			while( tid1 >= rd->nthreads && !(rd->error) && rd->nusers < rd->uids_uniq.size() );
			if(rd->error) break;
			__atomic_fetch_add(&(rd->running_threads),1,__ATOMIC_ACQ_REL);
			running = true;
		}
		
		unsigned int nusers = __atomic_fetch_add(&(rd->nusers),1,__ATOMIC_ACQ_REL);
		if(nusers >= rd->uids_uniq.size()) break;
		int64_t uid0 = rd->uids_uniq[nusers];
		
		//process the given user
		unsigned int ne = 0;
		for(crecords::iterator cit = rd->tr->get_uid_iterator_begin(uid0); cit != rd->tr->get_uid_iterator_end(uid0); ++cit) {
			const record& r = *cit;
			if( r.ts < rd->tmin || ( rd->tmax[rd->ntmax-1] > 0 && r.ts >= rd->tmax[rd->ntmax-1] ) ) continue; // leave out records outside the given time range
			
			uint32_t tmin1 = r.ts - rd->dt1;
			uint32_t tmax1 = r.ts + rd->dt1;
			if(r.startstop == 0) tmax1 = r.ts + rd->dt2; // start of a trip, use the transit-related time interval
			if(r.startstop == 1) tmin1 = r.ts - rd->dt2; // end of a trip, use the transit-related time interval
			if(tmin1 < rd->tmin) tmin1 = rd->tmin;
			if(rd->tmax[rd->ntmax-1] > 0 && tmax1 >= rd->tmax[rd->ntmax-1]) tmax1 = rd->tmax[rd->ntmax-1] - 1;
			
			if(rd->statall) {
				// search for any temporal matches
				if(rd->cr->search_records(r.cell_id, tmin1, tmax1, s_candidates, s_candidates, true)) {
					if(!rd->skiperrors) {
						fprintf(stderr,"Error searching for record with the new method: (%ld, %u, %d)!\n",r.uid,r.ts,r.cell_id);
						rd->error = true;
						break;
					}
				}
			}
			else {
				// search for only possible matches
				if(rd->cr->search_records_possible(r.cell_id, tmin1, tmax1, s_candidates)) {
					if(!rd->skiperrors) {
						fprintf(stderr,"Error searching for record with the new method: (%ld, %u, %d)!\n",r.uid,r.ts,r.cell_id);
						rd->error = true;
						return 0;
					}
				}
			}
		
			ne++;
		} //for(i) -- process all records by one user
		
		
		//test for impossible matches in the set, write out results
		if(ne > 0) {
			if(ne > nemax) {
				for(size_t i=0;i<nres;i++) {
					h_possible_matches[i].resize(ne,0); //note: new elements are added with the specified value (0)
					h_impossible_matches[i].resize(ne,0); //the clear() at the end ensures that the vectors at this point are empty, so all elements become zero
				}
				nemax = ne;
			}
			
			// everything is in s_candidates, need to check all
			// store results in the histograms
			for(auto it = s_candidates.cbegin(); it != s_candidates.cend(); ++it) {
				unsigned int ni = 0;
				int ni2 = 0;
				
				// check if user pair is among the histograms needed
				std::pair<uint32_t,uint32_t> nr;
				uint32_t ne2 = (uint32_t)(rd->cru->find_user(*it).second);
				if(rd->searchreverse) nr = std::make_pair(ne,ne2);
				else nr = std::make_pair(ne2,ne);
				std::vector<sbounds>::const_iterator it1 = std::lower_bound(vbounds.begin(),vbounds.end(),nr);
				if(it1 == vbounds.end()) continue;
				if( !((*it1).contains(nr.first,nr.second)) ) continue;
				size_t j1 = it1 - vbounds.begin();
				
				if(rd->ntmax > 1) {
					for(int i=0;i<rd->ntmax;i++) { matches_possible[i] = 0; matches_temporal[i] = 0; }
					if(rd->searchreverse) ni2 = rd->cru->compare_users_onematch_tmulti(*it,*(rd->tr),uid0,rd->ntmax,matches_temporal,matches_possible,
						rd->tmin,rd->tmax,false,true);
					else ni2 = rd->tr->compare_users_onematch_tmulti(uid0,*(rd->cru),*it,rd->ntmax,matches_temporal,matches_possible,
						rd->tmin,rd->tmax,false,true);
					if(ni2 == -1) {
						fprintf(stderr,"Error comparing users %ld and %ld!\n",uid0,*it);
						rd->error = true;
						break;
					}
					for(int i=0;i<rd->ntmax;i++) { // process result
						// add to the histograms
						if(matches_possible[i]) h_possible_matches[i*(vbounds.size()) + j1][ matches_possible[i]-1 ]++;
						if(matches_temporal[i] && rd->statall) h_impossible_matches[i*(vbounds.size()) + j1][ matches_temporal[i]-1 ]++;
					}
				}
				else {
					if(rd->searchreverse) {
						if(rd->onematch) ni2 = rd->cru->compare_users2_onematch(*it,*(rd->tr),uid0,&ni,0,false,rd->skiperrors,true,rd->tmin,rd->tmax[0]);
						else ni2 = rd->cru->compare_users2(*it,*(rd->tr),uid0,&ni,0,false,rd->skiperrors,true,rd->tmin,rd->tmax[0]);
					}
					else {
						if(rd->onematch) ni2 = rd->tr->compare_users2_onematch(uid0,*(rd->cru),*it,&ni,0,false,rd->skiperrors,true,rd->tmin,rd->tmax[0]);
						else ni2 = rd->tr->compare_users2(uid0,*(rd->cru),*it,&ni,0,false,rd->skiperrors,true,rd->tmin,rd->tmax[0]);
					}
					
					if(ni2 == -1) {
						fprintf(stderr,"Error comparing users %ld and %ld!\n",uid0,*it);
						rd->error = true;
						break;
					}
					
					// add to the histograms
					if(ni > 0) {
						if(rd->statall) h_impossible_matches[j1][ni-1]++;
						if(ni2 == 1) h_possible_matches[j1][ni-1]++;
					}
				}
			}
			
			s_candidates.clear();
		} //if(ne > 0)
		
		if(rd->out_progress_uid) if( nusers % (rd->out_progress_uid) == 0) { fprintf(stderr,"\r%u users processed",nusers); fflush(stderr); }
	
	} //while( ! error ) -- process all users
	
	// decrement the counter with the number of running threads
	tid1 = __atomic_fetch_sub(&(rd->running_threads),1,__ATOMIC_ACQ_REL);
	if(tid1 == 0) { // subtracting 1 from unsigned 0 is undefined behavior ?
		rd->running_threads = 0;
		rd->error = true;
	}
	
	// add up the main result histogram
	if( ! (rd->error) ) {
		pthread_mutex_lock(&(rd->res_mutex));
		for(size_t i=0;i<nres;i++) {
			if(rd->res_possible_matches[i].size() < h_possible_matches[i].size()) rd->res_possible_matches[i].resize(h_possible_matches[i].size(),0);
			if(rd->statall) if(rd->res_temporal_matches[i].size() < h_impossible_matches[i].size()) rd->res_temporal_matches[i].resize(h_impossible_matches[i].size(),0);
			
			for(size_t j=0;j<h_possible_matches[i].size();j++) rd->res_possible_matches[i][j] += h_possible_matches[i][j];
			if(rd->statall) for(size_t j=0;j<h_impossible_matches[i].size();j++) rd->res_temporal_matches[i][j] += h_impossible_matches[i][j];
		}
		pthread_mutex_unlock(&(rd->res_mutex));
	}
	
	delete[]h_possible_matches;
	delete[]h_impossible_matches;
	if(matches_possible) delete[]matches_possible;
	if(matches_temporal) delete[]matches_temporal;
	
	return 0;
}


// read the records and optionally filter the users in it
int read_records(crecords& r, const char* fn, bool inzip, const char* userfilter) {
	if(userfilter == 0) return crecords::read_crecords_serialized(r,fn,inzip);
	std::unordered_set<int64_t> users;
	FILE* uf = fopen(userfilter,"r");
	if(uf == 0) {
		fprintf(stderr,"read_records(): Error opening user list file %s!\n",userfilter);
		return 1;
	}
	while(1) {
		int64_t uid;
		int a = fscanf(uf,"%ld",&uid);
		if(a == EOF) break;
		if(a != 1) {
			fprintf(stderr,"read_records(): Error reading user ids from file %s!\n",userfilter);
			fclose(uf);
			return 1;
		}
		users.insert(uid);
	}
	fclose(uf);
	
	crecords tmp;
	int r0 = crecords::read_crecords_serialized(tmp,fn,inzip);
	if(r0) return r0;
	
	filter_iterator<crecords::const_iterator,crecords::sentinel> fi(tmp.begin(),tmp.end());
	fi.set_userfilter(&users);
	r.add_records_custom(fi,record_iterator_sentinel());
	return 0;
}

int main(int argc, char **argv)
{
	char** cdrfiles = 0; //files with CDR records (multiple binary files can be given, which should contain the same records with different sort orders)
	unsigned int ncdrfiles = 0;
	char* ltafile = 0; //file with transportation data
	char* ltauserlist = 0; // filter users in the lta data as given
	char* cdruserlist = 0; // filter users to the ones given here
	
	char** placesmapfiles = 0; //place-cell mapping files
	unsigned int nplacesmapfiles = 0; // number of files to read
	
	unsigned int dt1 = 600; // time window for walking (default: 10 min)
	unsigned int dt2 = 300; // time window for transit (default: 5 min)
	double radius1p = 500.0; // possible radius for walking (all of these need to be present among the place mappings given as placesmapfiles)
	double radius1i = 500.0; // impossible radius for walking
	double radius2p = 1000.0; // possible radius for transit
	double radius2i = 2000.0; // impossible radius for transit
		
	uint32_t tmin = 0; // minimum and maximum timestamp
	uint32_t* tmax = 0; // to limit searches in time
	uint32_t tmaxnull = 0; // to use as a limit if none is given
	int ntmax = 0; // several tmax values to get statistics of increasing intervals
	
	unsigned int nthreads = 1;
	
	bool searchreverse = false; //if set, search for CDR events among LTA records
	bool statall = false; //if true, all temporal matches are considered
	bool onematch = false; //if true, each point can result in only one match also in the searched dataset
	bool inzip = false; //read compressed binary files
	
	// auxilliary data for the custom popen() implementation used to save gziped files (note: the glibc / stdlib.h popen() will probably request
	//	to reserve memory 2x the size of the calling process due to using fork() which might fail)
	struct popen_noshell_pass_to_pclose spc,spc2;
	const char* popen_args[4];
	popen_args[0] = "/bin/gzip";
	popen_args[1] = "-c";
	popen_args[2] = 0;
	popen_args[3] = 0;
	
	bool tr_toverlap = true; // control wether time intervals can overlap when checking record matches
	bool cr_toverlap = true; //	(default: they can, but can be set to false if needed)
	
	std::vector<sbounds> vbounds; // bounds to limit search to
	bool bstdin = false; // read bounds from stdin
	char* bfile = 0; // read bounds parameters from this file
	
	unsigned int out_progress_uid = 0;
	
	//!! TODO: more rigorous processing of arguments, protect against invalid arguments
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'p':
			int j;
			for(j=i+1;j<argc;j++) if(argv[j][0] == '-') break;
			nplacesmapfiles = j-i-1;
			if(nplacesmapfiles > 0) placesmapfiles = argv+i+1;
			else fprintf(stderr,"Missing file names after %s!\n",argv[i]);
			break;
		case 'r':
			if(argv[i][2] == 'U') cdruserlist = argv[i+1];
			else {
				int j;
				for(j=i+1;j<argc;j++) if(argv[j][0] == '-') break;
				ncdrfiles = j-i-1;
				if(ncdrfiles > 0) cdrfiles = argv + i + 1; // save the actual filenames too
			}
			break;
		case 'l':
			if(argv[i][2] == 'U') ltauserlist = argv[i+1];
			else ltafile = argv[i+1];
			break;
		case 'z':
			inzip = true;
			break;
		case 't':
			if(argv[i][2] == 'i') { // time intervals for the search; expected as tmin, tmax1, tmax2, ...
				int j;
				for(j=i+1;j<argc;j++) if(argv[j][0] == '-') break;
				if(j-i-1 < 2) {
					fprintf(stderr,"Invalid parameter: %s (-ti requires at least two timestamps)!\n",argv[i]);
					break;
				}
				ntmax = j-i-2;
				tmax = (uint32_t*)malloc(sizeof(uint32_t)*ntmax);
				if(tmax == 0) {
					fprintf(stderr,"Error allocating memory!\n");
					return 1;
				}
				tmin = atoi(argv[i+1]);
				for(int k=i+2;k<j;k++) tmax[k-i-2] = atoi(argv[k]);
				break;
			}
			dt1 = atoi(argv[i+1]);
			if(i+2 < argc && argv[i+2][0] != '-') { dt2 = atoi(argv[i+2]); i++; }
			else dt2 = dt1;
			i++;
			break;
		case 'd':
			if(argv[i][2] == '1') {
				radius1p = atof(argv[i+1]);
				if(i+2 < argc && argv[i+2][0] != '-') { radius1i = atof(argv[i+2]); i++; }
				else radius1i = radius1p;
				i++;
				break;
			}
			if(argv[i][2] == '2') {
				radius2p = atof(argv[i+1]);
				if(i+2 < argc && argv[i+2][0] != '-') { radius2i = atof(argv[i+2]); i++; }
				else radius2i = radius2p;
				i++;
				break;
			}
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
		case 'T':
			nthreads = atoi(argv[i+1]);
			break;
		case 'R':
			searchreverse = true;
			break;
		case 'A':
			statall = true;
			break;
		case '1':
			onematch = true;
			break;
		case 'O':
			if(argv[i][2] == 'l') { tr_toverlap = false; break; }
			if(argv[i][2] == 'c') { cr_toverlap = false; break; }
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
		case 'b':
			// read all four bounds for generation at once -- or - means read these from the stdin (possible multiple variations)
			if(argv[i+1][0] == '-') { bstdin = true; i++; break; }
			else { bfile = argv[i+1]; i++; break; }
			break;
		case 'P':
			out_progress_uid = atoi(argv[i+1]);
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	
	if(ncdrfiles == 0 || cdrfiles == 0 || ltafile == 0 || placesmapfiles == 0) {
		fprintf(stderr,"Error: missing input files!\n");
		return 1;
	}
	
	rdata rd;
	crecords* cr1 = new crecords();
	crecords* cr2 = 0;
	if(ncdrfiles > 1) cr2 = new crecords();
	crecords tr;
	cr1->timeoverlap = cr_toverlap;
	if(cr2) cr2->timeoverlap = cr_toverlap;
	tr.timeoverlap = tr_toverlap;
	rd.out_progress_uid = out_progress_uid;
	
	time_t t1;
	char* ctimec = 0;
	// read bounds from stdin and sort them
	if(bstdin || bfile) { // read a list of bounds from stdin, store them in vectors and run for each
		FILE* bin = stdin;
		if(bfile) bin = fopen(bfile,"r");
		if(bin == 0) {
			fprintf(stderr,"Error opening file %s!\n",bfile);
			return 1;
		}
		while(1) {
			sbounds b;
			if(!b.readbounds(bin)) break;
			vbounds.push_back(b);
			int a;
			while(1) { a = getc(bin); if(a == '\n') break; if( ! (a == ' ' || a == '\t') ) break; }
			if(a != '\n') break;
		}
		if(vbounds.size() == 0) {
			fprintf(stderr,"Error: expected a list of bounds as input, none could be read!\n");
			return 1;
		}
		if(bfile) fclose(bin);
		t1 = time(0);
		ctimec = ctime(&t1);
		for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
		if(bfile) fprintf(stderr,"%s, %lu bounds parameters read from file %s\n",ctimec,vbounds.size(),bfile);
		else fprintf(stderr,"%s, %lu bounds parameters read from stdin\n",ctimec,vbounds.size());
	}
	if(vbounds.size() > 1) std::sort(vbounds.begin(),vbounds.end());
	
	//read place-cell mapping
	for(unsigned int i=0;i<nplacesmapfiles;i++) {
		FILE* pmf = fopen(placesmapfiles[i],"r");
		if(!pmf) {
			fprintf(stderr,"Error opening place-cell mapping file %s!\n",placesmapfiles[i]);
			return 1;
		}
		fprintf(stderr,"Reading place-cell mapping from file %s -- ",placesmapfiles[i]);
		if(cr1->read_places_map(pmf) != 0) {
			fprintf(stderr,"\nError reading place-cell mapping from file %s!\n",placesmapfiles[i]);
			fclose(pmf);
			return 1;
		}
		fclose(pmf);
	}
	tr.copy_places_map_reverse(*cr1); // copy everything reversed to the transportation records
	if(cr2) cr2->copy_places_map(*cr1); // copy into the second one as well
	
	// read CDR and LTA records first
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, reading CDR and transportation data\n",ctimec);
	for(unsigned int j=0;j<ncdrfiles && j<2;j++) {
		int r1;
		uint64_t r2;
		if(j==0) { r1 = read_records(*cr1,cdrfiles[j],inzip,cdruserlist); r2 = cr1->getnrecords(); }
		if(j==1) { r1 = read_records(*cr2,cdrfiles[j],inzip,cdruserlist); r2 = cr2->getnrecords(); }
		if(r1) {
			fprintf(stderr,"Error reading data from file %s!\n",cdrfiles[j]);
			return 2;
		}
		t1 = time(0);
		ctimec = ctime(&t1);
		for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
		fprintf(stderr,"%s, %lu CDR records read from input file %s\n",ctimec,r2,cdrfiles[j]);
	}
	if(read_records(tr,ltafile,inzip,ltauserlist)) {
		fprintf(stderr,"Error reading data from file %s!\n",ltafile);
		return 2;
	}
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, %lu transportation records read from input file %s\n",ctimec,tr.getnrecords(),ltafile);
	uint64_t lines2 = 0;
	
	
	// create copies
	crecords* cru = 0;
	double radiusp = radius1p;
	if(radius2p > radius1p) radiusp = radius2p;
	if(searchreverse) {
		cru = new crecords(tr);
		cru->set_compare_params(dt1,dt2,radius1p,radius1i,radius2p,radius2i); // create place -- cell mappings for user comparisons
		tr.set_default_places_map(radiusp);
		
		rd.cr = &tr;
		if( (cr1->get_sortorder() == crecords::sortorder::uid_time || cr1->get_sortorder() == crecords::sortorder::uid_time_cellid) || ncdrfiles == 1) {
			rd.tr = cr1;
			delete cr2;
			cr2 = 0;
		}
		else {
			rd.tr = cr2;
			delete cr1;
			cr1 = 0;
		}
		rd.cru = cru;
		rd.searchreverse = true;
		lines2 = cr1->getnrecords();
	}
	else {
		// cr1: should be ordered by cell_id, time
		// cr2: should be ordered by user_id, time
		if(cr2 == 0) cr2 = new crecords(*cr1); // just copy everything
		else {
			if( cr1->get_sortorder() == crecords::sortorder::uid_time || cr1->get_sortorder() == crecords::sortorder::uid_time_cellid ||
				cr2->get_sortorder() == crecords::sortorder::cellid_time || cr2->get_sortorder() == crecords::sortorder::cellid_time_uid ) {
					// swap the two pointers
					crecords* tmp = cr1;
					cr1 = cr2;
					cr2 = tmp;
				}
		}
		
		tr.set_compare_params(dt1,dt2,radius1p,radius1i,radius2p,radius2i,true); // create place -- cell mappings for user comparisons
			// note: last param: places_map is supposed to be in "reverse", i.e. cell -> place as it is needed in this case
		cr1->set_default_places_map(radiusp);
		rd.cr = cr1;
		rd.cru = cr2;
		rd.tr = &tr;
		rd.searchreverse = false;
		lines2 = tr.getnrecords();
	}
	
	rd.dt1 = dt1;
	rd.dt2 = dt2;
	rd.skiperrors = true;
	rd.statall = statall;
	rd.onematch = onematch;
	//~ rd.statall_possible = statall_possible;
	rd.tmin = tmin;
	if(ntmax > 0) {
		rd.tmax = tmax;
		rd.ntmax = ntmax;
	}
	else {
		rd.tmax = &tmaxnull;
		rd.ntmax = 1;
	}
	rd.error = false;
	
	
	
	//explicit sort step is needed -- note: if data is already sorted, this is a no-op
	//do this on separate threads per sort step
	//note: nested functions is probably a gcc extension
	{	
		struct pthread_sort {
			static void* pthread_sort_by_cellid_time(void* p) { crecords* c = (crecords*)p; c->sort_by_cellid_time(); return 0; }
			static void* pthread_sort_by_uid_time(void* p) { crecords* c = (crecords*)p; c->sort_by_uid_time(); return 0; }
			static void* pthread_sort_by_time(void* p) { crecords* c = (crecords*)p; c->sort_by_time(); return 0; }
		};
		pthread_t sort_threads[3];
		bool running_threads[3] = { false, false, false };
		
		if(statall) {
			if(rd.cr->get_sortorder() != crecords::sortorder::time) {
				if(nthreads > 1) {
					if(pthread_create(sort_threads,0,&pthread_sort::pthread_sort_by_time,(void*)(rd.cr) )) { fprintf(stderr,"Error creating sort threads!\n"); return 6; }
					running_threads[0] = true;
				}
				else rd.cr->sort_by_time();
			}
		}
		else {
			if( ! (rd.cr->get_sortorder() == crecords::sortorder::cellid_time || rd.cr->get_sortorder() == crecords::sortorder::cellid_time_uid) ) {
				if(nthreads > 1) {
					if(pthread_create(sort_threads,0,&pthread_sort::pthread_sort_by_cellid_time,(void*)(rd.cr) )) { fprintf(stderr,"Error creating sort threads!\n"); return 6; }
					running_threads[0] = true;
				}
				else rd.cr->sort_by_cellid_time();
			}
			else rd.cr->create_cell_idx();
		}
		if( ! (rd.tr->get_sortorder() == crecords::sortorder::uid_time || rd.tr->get_sortorder() == crecords::sortorder::uid_time_cellid) ) {
			if(nthreads > 1) {
				if(pthread_create(sort_threads + 1,0,&pthread_sort::pthread_sort_by_uid_time,(void*)(rd.tr) )) { fprintf(stderr,"Error creating sort threads!\n"); return 6; }
				running_threads[1] = true;
			}
			else rd.tr->sort_by_uid_time();
		}
		else rd.tr->create_uid_idx();
		if( ! (rd.cru->get_sortorder() == crecords::sortorder::uid_time || rd.cru->get_sortorder() == crecords::sortorder::uid_time_cellid) ) {
			if(nthreads > 1) {
				if(pthread_create(sort_threads + 2,0,&pthread_sort::pthread_sort_by_uid_time,(void*)(rd.cru) )) { fprintf(stderr,"Error creating sort threads!\n"); return 6; }
				running_threads[2] = true;
			}
			else rd.cru->sort_by_uid_time();
		}
		else rd.cru->create_uid_idx();
		
		if(running_threads[0]) if(pthread_join(sort_threads[0],0)) { fprintf(stderr,"Error runnning sort threads!\n"); return 7; }
		if(running_threads[1]) if(pthread_join(sort_threads[1],0)) { fprintf(stderr,"Error runnning sort threads!\n"); return 7; }
		if(running_threads[2]) if(pthread_join(sort_threads[2],0)) { fprintf(stderr,"Error runnning sort threads!\n"); return 7; }
	}
	
	//fill out start offsets for all users
	int64_t uid0 = 0;
	int64_t uid1 = 0;
	
	rd.uids_uniq = rd.tr->get_uids(false);
	
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, done preprocessing, starting search\n",ctimec);
	
	// allocate results vectors
	const size_t nres = rd.ntmax * vbounds.size();
	rd.res_possible_matches = new std::vector<int64_t>[nres];
	rd.res_temporal_matches = new std::vector<int64_t>[nres];
	rd.vbounds = &vbounds;
	pthread_mutex_init(&(rd.res_mutex),0);
	
	rd.nusers = 0;
	rd.running_threads = nthreads;
	rd.nthreads = nthreads;
	rd.tid = 0;
	if(nthreads > 1) {
		pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t)*nthreads);
		if(!threads) {
			fprintf(stderr,"Error allocating memory!\n");
			return 4;
		}
	
		//create threads, wait until they run
		unsigned int i=0;
		for(i=0;i<nthreads;i++) {
			if(pthread_create(threads + i,0,&process_records,(void*)&rd)) {
				fprintf(stderr,"Error creating thread %u!\n",i);
				rd.error = true;
			}
		}
		rd.running_threads = i;
		rd.nthreads = i;
		
		// wait until the threads are finished or SIGINT is received in which case the user is asked to adjust the number of running threads
		sigset_t set;
		sigemptyset(&set);
		sigaddset(&set, SIGQUIT);
		sigaddset(&set, SIGINT);
		pthread_sigmask(SIG_BLOCK,&set,0);
		while( rd.running_threads && !rd.error ) {
			struct timespec ts;
			ts.tv_sec = 5;
			ts.tv_nsec = 0;
			errno = 0;
			int r = sigtimedwait(&set,0,&ts);
			if(r > 0) {
				pthread_sigmask(SIG_UNBLOCK,&set,0);
				fprintf(stderr,"Number of running threads: %u, possible range: 1 -- %u\nSpecify the new value in this range or 0 to quit: ",
					rd.running_threads,i);
				fflush(stderr);
				unsigned int nt2;
				r = fscanf(stdin,"%u",&nt2);
				if(r != 1) fprintf(stderr,"No number given, not adjusting\n");
				else {
					if(nt2 > i) fprintf(stderr,"Error: out of range, not adjusting\n");
					else {
						if(nt2 == 0) fprintf(stderr,"Exiting...\n");
						else fprintf(stderr,"Limiting the number of running threads to %u\n",nt2);
						rd.nthreads = nt2;
					}
				}
				pthread_sigmask(SIG_BLOCK,&set,0);
			}
			else if(errno != EAGAIN) {
				fprintf(stderr,"Error calling sigtimedwait(), aborting...\n");
				rd.error = true;
				break;
			}
		}
		pthread_sigmask(SIG_UNBLOCK,&set,0);
		
		for(unsigned int j=0;j<i;j++) {
			if(pthread_join(threads[j],0)) {
				fprintf(stderr,"Error waiting for thread %u!\n",i);
			}
		}
		
		free(threads);
	}
	else process_records((void*)&rd); // run everything in the current thread if we are not using > 1 threads to make debugging simpler
	
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	unsigned int nusers1 = rd.nusers;
	if(nusers1 > rd.uids_uniq.size()) nusers1 = rd.uids_uniq.size();
	fprintf(stderr,"\n%s, %u users processed\n",ctimec,nusers1);
	
	// write out results (everything to stdout)
	FILE* fout = stdout;
	for(int i=0;i<ntmax;i++) for(size_t j=0;j<vbounds.size();j++) {
		size_t k = vbounds.size()*i + j;
		size_t s1 = rd.res_possible_matches[k].size();
		if(rd.res_temporal_matches[k].size() > s1) s1 = rd.res_temporal_matches[k].size();
		for(size_t l=0;l<s1;l++) {
			vbounds[j].printbounds(fout);
			fprintf(fout,"\t%u\t%u\t%lu",tmin,tmax[i],l+1);
			if(l < rd.res_possible_matches[k].size()) fprintf(fout,"\t%lu",rd.res_possible_matches[k][l]);
			else fprintf(fout,"\t0");
			if(l < rd.res_temporal_matches[k].size()) fprintf(fout,"\t%lu",rd.res_temporal_matches[k][l]);
			else fprintf(fout,"\t0");
			fputc('\n',fout);
		}
	}
	
	delete[]rd.res_possible_matches;
	delete[]rd.res_temporal_matches;
	pthread_mutex_destroy(&rd.res_mutex);
	
	if(cru) delete cru;
	if(cr1) delete cr1;
	if(cr2) delete cr2;
	free(tmax);
	
	return 0;
}

