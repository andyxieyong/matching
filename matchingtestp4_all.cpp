/*
 * matchingtestp4_all.cpp
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
 * 	count impossible matches (i.e.~total number of temporal matches), write out statistics about that too
 * 	keep track of all points, allow only one match per point in both datasets
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
#include <algorithm> // std::sort for determining the top matches
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include "cdrrecords.h"
#include "popen_noshell.h"


// shared data by the threads (records to process) -- note: this could be implemented as a derived class of crecords
struct rdata {
	//distinct user ids and pointers in the above vectors to the start of the data associated with them
	std::vector<int64_t> uids_uniq;
	//~ std::vector<uint64_t> uids_off;
	
	//number of users processed (a thread should choose it's next user to process as uids_uniq[nusers]
	unsigned int nusers;
	//~ pthread_mutex_t nusers_mutex; //mutex to use to choose the next user to process
	
	//files to write the output to (possible matches and summary statistics in separate files)
	FILE* out_possible;
	FILE* out_summary;
	
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
	
	// 4. optional: store the CDR data again, sorted by timestamp to do the original kind of processing to check if the new version is working well
	//~ crecords* crc; -- note: removed these so as to simplify the process_records() function
	
	//if set to true, failed tests do not result in aborted run, so that these can be better inspected from the debugger
	//~ bool debugrun;
	
	bool skiperrors; //if true, errors from search_records_possible() (i.e. missing place ID) do not result in an aborted run
		//when searching in reverse (searching for CDR records among the LTA data), there can be cells with no stops in or around them,
		//that is not an error
	bool statall; //if true, result (output) is a statistic of the number of all temporal matches, separated to possible and impossible
	bool statall_possible; //if true, the summary file will contain only the statistic of possible matches (i.e. all temporal matches need not be calculated)
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
	std::unordered_set<int64_t> s_candidates; //use this to store all / any matches which are later checked 
	
	std::vector<uint32_t>* h_possible_matches = new std::vector<uint32_t>[rd->ntmax]; //used as temporary storage to create histograms of the number of matches
	std::vector<uint32_t>* h_impossible_matches = new std::vector<uint32_t>[rd->ntmax]; //note: these are resized dynamically as needed
	unsigned int* matches_possible = 0;
	unsigned int* matches_temporal = 0;
	if(rd->ntmax > 1) {
		matches_possible = new unsigned int[rd->ntmax];
		matches_temporal = new unsigned int[rd->ntmax];
	}
	
	// assign thread ID
	unsigned int tid1 = __atomic_fetch_add(&(rd->tid),1,__ATOMIC_ACQ_REL);
	bool running = true;
	
	while(1) {
		// check for error or quit request
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
			
			if(rd->cr->search_records(r.cell_id, tmin1, tmax1, s_candidates, s_candidates, true)) {
				if(!rd->skiperrors) {
					fprintf(stderr,"Error searching for record with the new method: (%ld, %u, %d)!\n",r.uid,r.ts,r.cell_id);
					rd->error = true;
					break;
				}
			}
			ne++;
		} //for(i) -- process all records by one user
		if(rd->error) break;
		
		//test for impossible matches in the set, write out results
		if(ne > 0) {
			for(int i=0;i<rd->ntmax;i++) {
				h_possible_matches[i].resize(ne,0); //note: new elements are added with the specified value (0)
				h_impossible_matches[i].resize(ne,0); //the clear() at the end ensures that the vectors at this point are empty, so all elements become zero
			}
			
			for(auto it = s_candidates.cbegin(); it != s_candidates.cend(); ++it) {
				// check for possible matches first
				unsigned int ni = 0;
				int ni2 = 0;
				
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
					for(int i=0;i<rd->ntmax;i++) {
						if(matches_possible[i]) h_possible_matches[i][ matches_possible[i]-1 ]++;
						if(matches_temporal[i]) h_impossible_matches[i][ matches_temporal[i]-1 ]++;
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
					
					if(ni2 == 0) { // impossible, ni contains the number of temporal matches
						// add to the histogram
						if(ni > 0) h_impossible_matches[0][ni-1]++;
					}
					else { // possible
						// just add to the histogram
						if(ni > 0) h_possible_matches[0][ni-1]++;
					}
				}
			}
			
			// output histograms (note: do not output matches, that can be done with the original matchingtest4.cpp program which does the search
			//			for only possible matches more efficiently)
			flockfile(rd->out_summary);
			
			for(uint32_t i=0;i<ne;i++) {
				if(rd->ntmax > 1) {
					for(int j=0;j<rd->ntmax;j++) if(h_possible_matches[j][i] || h_impossible_matches[j][i]) 
						fprintf(rd->out_summary,"%u\t%u\t%ld\t%u\t%u\t%u\t%u\n",rd->tmin,rd->tmax[j],uid0,ne,i+1,
							h_possible_matches[j][i],h_impossible_matches[j][i]);
				}
				else if(h_possible_matches[0][i] || h_impossible_matches[0][i])
					fprintf(rd->out_summary,"%ld\t%u\t%u\t%u\t%u\n",uid0,ne,i+1,h_possible_matches[0][i],h_impossible_matches[0][i]);
			}
			
			funlockfile(rd->out_summary);
			
			for(int i=0;i<rd->ntmax;i++) {
				h_possible_matches[i].clear(); //this ensures that all counts will be zero the next time
				h_impossible_matches[i].clear(); //this ensures that all counts will be zero the next time
			}
				
			s_candidates.clear();
				
		} //if(ne > 0)
	
	} //while( ! error ) -- process all users
	
	// decrement the counter with the number of running threads
	tid1 = __atomic_fetch_sub(&(rd->running_threads),1,__ATOMIC_ACQ_REL);
	if(tid1 == 0) { // subtracting 1 from unsigned 0 is undefined behavior ?
		rd->running_threads = 0;
		rd->error = true;
	}
	
	delete[]h_possible_matches;
	delete[]h_impossible_matches;
	if(matches_possible) delete[]matches_possible;
	if(matches_temporal) delete[]matches_temporal;
	
	return 0;
}


int main(int argc, char **argv)
{
	char* cdrfile = 0; //file with CDR records
	char* ltafile = 0; //file with transportation data
	
	char** placesmapfiles = 0; //place-cell mapping files
	unsigned int nplacesmapfiles = 0; // number of files to read
	
	unsigned int dt1 = 600; // time window for walking (default: 10 min)
	unsigned int dt2 = 300; // time window for transit (default: 5 min)
	double radius1p = 500.0; // possible radius for walking
	double radius1i = 500.0; // impossible radius for walking
	double radius2p = 1000.0; // possible radius for transit
	double radius2i = 2000.0; // impossible radius for transit
	
	char* outf_summary = 0; //write out the number of impossible matches and excluded users for each input user separately instead of in the main output file
	
	uint32_t tmin = 0; // minimum and maximum timestamp
	uint32_t* tmax = 0; // to limit searches in time
	uint32_t tmaxnull = 0; // to use as a limit if none is given
	int ntmax = 0; // several tmax values to get statistics of increasing intervals
	
	//Voronoi tessellation
	char* v_vertex_coords = 0;
	char* v_vertex_to_cell = 0;
	char* v_edges = 0;
	char* v_cell_unique = 0;
	
	unsigned int nthreads = 1;
	
	bool searchreverse = false; //if set, search for CDR events among LTA records
	bool onematch = false; //if true, each point can result in only one match also in the searched dataset
	bool outzip = false; //compress output files (if not writing to stdout)
	bool inzip = false; //read compressed binary files
	
	// auxilliary data for the custom popen() implementation used to save gziped files (note: the glibc / stdlib.h popen() will probably request
	//	to reserve memory 2x the size of the calling process due to using fork() which might fail)
	struct popen_noshell_pass_to_pclose spc,spc2;
	const char* popen_args[4];
	popen_args[0] = "/bin/gzip";
	popen_args[1] = "-c";
	popen_args[2] = 0;
	popen_args[3] = 0;
	
	
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
			cdrfile = argv[i+1];
			break;
		case 'l':
			ltafile = argv[i+1];
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
		case '1':
			onematch = true;
			break;
		case 'o':
			if(argv[i][2] == 'z') { outzip = true; break; }
			outf_summary = argv[i+1];
			break;
		case 'v':
			if(argv[i][2] == 'c') { v_vertex_coords = argv[i+1]; break; }
			if(argv[i][2] == 'v') { v_vertex_to_cell = argv[i+1]; break; }
			if(argv[i][2] == 'e') { v_edges = argv[i+1]; break; }
			if(argv[i][2] == 'u') { v_cell_unique = argv[i+1]; break; }
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	
	if(cdrfile == 0 || ltafile == 0 || placesmapfiles == 0) {
		fprintf(stderr,"Error: missing input files!\n");
		return 1;
	}
	
	if(ntmax > 1 && onematch == false) {
		fprintf(stderr,"Error: searching with multiple time windows without limiting to one match per point is not supported (use the -1 switch)!\n");
		free(tmax);
		return 1;
	}
		
	rdata rd;
	crecords cr;
	crecords tr;
	cr.use_combined_sort = true;
	tr.use_combined_sort = true;
	if(searchreverse) {
		rd.cr = &tr;
		rd.tr = &cr;
		rd.searchreverse = true;
	}
	else {
		rd.cr = &cr;
		rd.tr = &tr;
		rd.searchreverse = false;
	}
	
	rd.dt1 = dt1;
	rd.dt2 = dt2;
	rd.skiperrors = true;
	rd.onematch = onematch;
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
	
	//read place-cell mapping
	for(unsigned int i=0;i<nplacesmapfiles;i++) {
		FILE* pmf = fopen(placesmapfiles[i],"r");
		if(!pmf) {
			fprintf(stderr,"Error opening place-cell mapping file %s!\n",placesmapfiles[i]);
			free(tmax);
			return 1;
		}
		fprintf(stderr,"Reading place-cell mapping from file %s -- ",placesmapfiles[i]);
		if(cr.read_places_map(pmf) != 0) {
			fprintf(stderr,"\nError reading place-cell mapping from file %s!\n",placesmapfiles[i]);
			free(tmax);
			fclose(pmf);
			return 1;
		}
		fclose(pmf);
	}
	tr.copy_places_map_reverse(cr); // copy everything reversed to the transportation records

	// read CDR and LTA records first
	time_t t1 = time(0);
	char* ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, reading CDR and transportation data\n",ctimec);
	if(crecords::read_crecords_serialized(cr,cdrfile,inzip)) {
		fprintf(stderr,"Error reading data from file %s!\n",cdrfile);
		free(tmax);
		return 2;
	}
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, %lu CDR records read from input file %s\n",ctimec,cr.getnrecords(),cdrfile);
	if(crecords::read_crecords_serialized(tr,ltafile,inzip)) {
		fprintf(stderr,"Error reading data from file %s!\n",ltafile);
		free(tmax);
		return 2;
	}
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, %lu transportation records read from input file %s\n",ctimec,tr.getnrecords(),ltafile);
	uint64_t lines2 = 0;
	if(searchreverse) lines2 = cr.getnrecords();
	else lines2 = tr.getnrecords();
	
	
	// open output files
	if(outf_summary) {
		if(outzip) {
			popen_args[1] = "-c";
			popen_args[2] = 0;
			rd.out_summary = popen_noshell_redirect(popen_args[0],popen_args,"w",&spc,0,outf_summary,0);
		}
		else rd.out_summary = fopen(outf_summary,"w");
		if(!rd.out_summary) {
			fprintf(stderr,"Error opening output file %s!\n",outf_summary);
			free(tmax);
			return 1;
		}
	}
	else rd.out_summary = stdout;
	
	
	// create copies
	crecords* cru = 0;
	double radiusp = radius1p;
	if(radius2p > radius1p) radiusp = radius2p;
	if(searchreverse) {
		cru = new crecords(tr);
		cru->set_compare_params(dt1,dt2,radius1p,radius1i,radius2p,radius2i); // create place -- cell mappings for user comparisons
		tr.set_default_places_map(radiusp);
	}
	else {
		cru = new crecords(cr);
		tr.set_compare_params(dt1,dt2,radius1p,radius1i,radius2p,radius2i,true); // create place -- cell mappings for user comparisons
			// note: last param: places_map is created in "reverse", i.e. cell -> place as it is needed in this case
		cr.set_default_places_map(radiusp);
	}
	cru->use_combined_sort = true;
	rd.cru = cru;
	
	//explicit sort step is needed
	//do this on separate threads per sort step
	//note: nested functions is probably a gcc extension
	if(nthreads > 1) {
		struct pthread_sort {
			static void* pthread_sort_by_cellid_time(void* p) { crecords* c = (crecords*)p; c->sort_by_cellid_time(); return 0; }
			static void* pthread_sort_by_uid_time(void* p) { crecords* c = (crecords*)p; c->sort_by_uid_time(); return 0; }
			static void* pthread_sort_by_time(void* p) { crecords* c = (crecords*)p; c->sort_by_time(); return 0; }
		};
		
		pthread_t sort_threads[3];
		
		if(pthread_create(sort_threads,0,&pthread_sort::pthread_sort_by_time,(void*)(rd.cr) )) { fprintf(stderr,"Error creating sort threads!\n"); return 6; }
		if(pthread_create(sort_threads+1,0,&pthread_sort::pthread_sort_by_uid_time,(void*)(rd.tr) )) { fprintf(stderr,"Error creating sort threads!\n"); return 6; }
		if(pthread_create(sort_threads+2,0,&pthread_sort::pthread_sort_by_uid_time,(void*)(rd.cru) )) { fprintf(stderr,"Error creating sort threads!\n"); return 6; }
		
		if(pthread_join(sort_threads[0],0)) { fprintf(stderr,"Error runnning sort threads!\n"); return 7; }
		if(pthread_join(sort_threads[1],0)) { fprintf(stderr,"Error runnning sort threads!\n"); return 7; }
		if(pthread_join(sort_threads[2],0)) { fprintf(stderr,"Error runnning sort threads!\n"); return 7; }
	}
	else {
		//do sort only using the current thread
		rd.cr->sort_by_time();
		rd.tr->sort_by_uid_time();
		rd.cru->sort_by_uid_time();
	}
	
	//fill out start offsets for all users
	int64_t uid0 = 0;
	int64_t uid1 = 0;
	
	crecords* crutmp = rd.tr;
	rd.uids_uniq = crutmp->get_uids(false);
	
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, done preprocessing, starting search\n",ctimec);
	
	
	rd.nusers = 0;
	rd.running_threads = nthreads;
	rd.nthreads = nthreads;
	rd.tid = 0;
	if(nthreads > 1) {
		pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t)*nthreads);
		if(!threads) {
			fprintf(stderr,"Error allocating memory!\n");
			free(tmax);
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
		
		// wait until the threads are finished or SIGINT is received in which case the user is asked to lower the number of running threads
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
	fprintf(stderr,"%s, %u users processed\n",ctimec,nusers1);
	if(rd.out_summary != stdout) {
		if(outzip) pclose_noshell(&spc);
		else fclose(rd.out_summary);
	}
	
	delete cru;
	free(tmax);
	
	return 0;
}

