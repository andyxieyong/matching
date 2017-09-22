/*
 * sertest.cpp -- test saving / loading crecords classes in binary format (using the cereal library)
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


#include "cdrrecords.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fstream>
#include <unordered_set>
#include <cereal/archives/portable_binary.hpp>

#include "popen_noshell.h"


int main(int argc, char **argv)
{
	char** cdrfiles = 0; //files with CDR records
	unsigned int ncdrfiles = 0;
	bool cdrfileszip = false;
	unsigned int cdr_header_skip = 1; //header in CDR record files

	char* sfnin = 0; //file name for serialized data for loading
	char* sfnout = 0; //file name for serialized data for saving
	bool sinzip = false; //read gzipped serialized data (output is never zipped (due to lazyness), it has to be compressed separately later
	
	bool randomtest = false; //create completely random data -- note: generated IDs and timestamps are completely random, so this makes little sense;
	unsigned int seed = 1;		//	can be mainly used to test the speed of the read, write and sort operations
	uint64_t rnum = 100000000; //test on this many random numbers
	bool read_startstop = false; // read start / stop info for transportation records
	
	/* program operation (with respect to inputs):
	  at least one text file (cdrfiles, with -r) or a binary file (sfnin, with -l) or the random flag (-R) has to be specified; at most two of these can be present
	  1. both of these are read (text file -> cr; binary file -> cr2); if -R is given instead of one of these, random data is generated there
	  2. optionally a user list and / or cell list can be specified, and the inputs are limited to these users (-U flag)
	  3. if text files were specified and Voronoi-tesselation and place coords was also given, the place-cell mapping is constructed (i.e. in cr)
	  4. if text files were NOT specified, but a Voronoi-tesselation and place coords were given, the place-cell mapping in the loaded binary
			data (i.e. in cr2) is discarded and a new mapping is created based on these parameters
	  5. input (both cr and cr2 if given) is sorted according to any sort parameter specified
	  6. if both text and binary input was given, these are compared, and the program aborts with error if these differ
	  7. if output file name was specified (sfnout, with -s), data is saved in binary format at the given file
	  
	  note: from the two inputs and the output (-r, -l and -s parameters) at least two need to be given to make sense running this program
		(if only one input and no output is given, the program can still be used to measure how long it takes to read and sort data)
	  
	  typical operating modes:
	  1. text file -> binary conversion (with optional sort step): -r and -s need to be specified
	  2. compare text file and binary file: -r and -l need to be specified
	  3. change sort order and / or place -- cell mapping in binary file: -l and -s need to be specified
		in all cases the Voronoi-tesselation and places (-v*, -c and -p) and the sort flags (-S*) can be optionally specified
	
	  notes about the place -- cell matching:
	    if only the Voronoi-tesselation is given (-v* and -c but no -p), the cells are matched among themselves
	    if only places are given (-p but no -v* and -c), the places are matched among themselves
	    if both all given (-p, -c and -v*), the matching between places and cells is created
	    the type is saved in the crecords for newly created serialized data
	*/
	
	bool sort_by_uid_time = false; //sort the input by uid and time
	bool sort_by_time = false; //sort the input by time
	bool sort_by_cellid_time = false; //sort the input by cellid and time
	bool sort_by_cellid_time_uid = false; //sort by cellid and time and also uid (so that the sorted result should be unique and not depend on the sort algorithm used)
	
	char* ctimec = 0;
	time_t t1;
	
	bool cmp_ignore_places = false; // ignore places_maps and related data when comparing
	
	uint64_t nrecords = 0; // if given, reserve this much space for the records before reading input (to avoid trying to allocate too much memory
		//	due to std::vector doubling storage on autogrowth)
		
	char* userlist = 0; // file to read user list to write only a limited set of users
	char* celllist = 0; // file to read cell list from to limit cells
	uint32_t tmin = 0; // timestamps to limit records to
	uint32_t tmax = 0;
	
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'r':
			if(i+1 == argc) {
				ncdrfiles = 1; // just signal that this flag was given, data will be read from stdin
			}
			else {
				int j;
				for(j=i+1;j<argc;j++) if(argv[j][0] == '-') break;
				ncdrfiles = j-i-1;
				if(ncdrfiles == 0) ncdrfiles = 1; // just signal that this flag was given, data will be read from stdin
				else cdrfiles = argv + i + 1; // save the actual filenames too
			}
			break;
		case 'z':
			cdrfileszip = true;
			break;
		case 'Z':
			sinzip = true;
			break;
		case 'l':
			sfnin = argv[i+1];
			i++;
			break;
		case 's':
			sfnout = argv[i+1];
			i++;
			break;
		case 'N':
			nrecords = strtoull(argv[i+1],0,10);
			break;
		case 'R':
			randomtest = true;
			if(i+1 < argc && argv[i+1][0] != '-') rnum = strtoul(argv[i+1],0,10);
			break;
		case 'U':
			userlist = argv[i+1];
			i++;
			break;
		case 'C':
			celllist = argv[i+1];
			i++;
			break;
		case 'T':
			if(i+2 >= argc || ! (isdigit(argv[i+1][0]) && isdigit(argv[i+2][0]) )  ) { fprintf(stderr,"Invalid parameter: -T requires two positive numbers!\n"); break; }
			tmin = atoi(argv[i+1]);
			tmax = atoi(argv[i+2]);
			if(tmax <= tmin) fprintf(stderr,"Invalid parameters: %s %s %s!\n",argv[i],argv[i+1],argv[i+2]);
			else i += 2;
			break;
		case 't':
			read_startstop = true;
			break;
		case 'S':
			if(argv[i][2] == 'u') { sort_by_uid_time = true; break; }
			if(argv[i][2] == 't') { sort_by_time = true; break; }
			if(argv[i][2] == 'c') { sort_by_cellid_time = true; break; }
			if(argv[i][2] == 'C') { sort_by_cellid_time_uid = true; break; }
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	
	if(sfnin == 0 && sfnout == 0 && ncdrfiles == 0 && randomtest == false) {
		fprintf(stderr,"No input or output specified!\n");
		return 1;
	}
	if(sfnin && ncdrfiles && randomtest) {
		fprintf(stderr,"Error: reading text and binary and also generating random data is not supported!\n");
		return 1;
	}
	
	// read user list if needed -- the file should only contain user ids
	std::unordered_set<int64_t> userset;
	std::unordered_set<int> cellset;
	if(userlist) {
		FILE* fu = fopen(userlist,"r");
		if(!fu) {
			fprintf(stderr,"Error opening file %s!\n",userlist);
			return 1;
		}
		while( ! (feof(fu) || ferror(fu) ) ) {
			int64_t u;
			int a = fscanf(fu,"%ld",&u);
			if(a == EOF) break;
			if(a != 1) {
				fprintf(stderr,"Invalid user id in file %s!\n",userlist);
				fclose(fu);
				return 2;
			}
			userset.insert(u);
		}
		fclose(fu);
	}
	if(celllist) {
		FILE* fu = fopen(celllist,"r");
		if(!fu) {
			fprintf(stderr,"Error opening file %s!\n",celllist);
			return 1;
		}
		while( ! (feof(fu) || ferror(fu) ) ) {
			int u;
			int a = fscanf(fu,"%d",&u);
			if(a == EOF) break;
			if(a != 1) {
				fprintf(stderr,"Invalid user id in file %s!\n",celllist);
				fclose(fu);
				return 2;
			}
			cellset.insert(u);
		}
		fclose(fu);
	}
	
	
	// define crecords variables for later use
	crecords* cr = new crecords(); // data loaded from text file (if given)
	crecords* cr2 = new crecords(); // data loaded from binary file (if given)
	crecords* crrandom = 0; // data to generate random records in (either one of the previous)
	crecords* crsave = 0; // data to save records from (either one of the previous)
	
	// 1. load data from text files (if specified)
	if(ncdrfiles) {
		if(nrecords > 0) cr->reserve(nrecords);
		
		// 1.2 CDR / transportation data records
		for(unsigned int i = 0;i < ncdrfiles;i++) {
			FILE* f = 0;
			struct popen_noshell_pass_to_pclose spc;
			if(cdrfiles) {
				const char* const popen_args[] = { "/bin/gzip", "-cd", cdrfiles[i], 0 };
				
				if(cdrfileszip) {
					f = popen_noshell(popen_args[0],popen_args,"r",&spc,0);
				}
				else f = fopen(cdrfiles[i],"r");
				if(!f) {
					fprintf(stderr,"Could not open input file %s!\n",cdrfiles[i]);
					return 5;
				}
			}
			else {
				f = stdin;
				cdr_header_skip = 0;
			}
			
			uint64_t lines = 0;
			tsv_iterator ti(f,cdr_header_skip,0,read_startstop);
			if(userlist || celllist || (tmin > 0 && tmax > tmin) ) {
				filter_iterator<tsv_iterator,record_iterator_sentinel> fi(ti,record_iterator_sentinel());
				if(userlist) fi.set_userfilter(&userset);
				if(celllist) fi.set_cellsfilter(&cellset);
				if(tmin > 0 && tmax > tmin) fi.set_tsfilter(tmin,tmax);
				
				lines = cr->add_records_custom(fi,record_iterator_sentinel(),0,false);
			}
			else lines = cr->add_records_custom(ti,record_iterator_sentinel(),0,false);
			
			t1 = time(0);
			ctimec = ctime(&t1);
			for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
			if(cdrfiles) {
				fprintf(stderr,"%s, %lu CDR records read from input file %s\n",ctimec,lines,cdrfiles[i]);
				if(cdrfileszip) pclose_noshell(&spc);
				else fclose(f);
			}
			else fprintf(stderr,"%s, %lu CDR records read from standard input\n",ctimec,lines);
		}
		
		// 1.4 sort the input if needed
		if(sort_by_time || sort_by_uid_time || sort_by_cellid_time || sort_by_cellid_time_uid) {
			fprintf(stderr,"\tsorting records\n");
			if(sort_by_time) cr->sort_by_time();
			if(sort_by_uid_time) cr->sort_by_uid_time();
			if(sort_by_cellid_time) cr->sort_by_cellid_time();
			if(sort_by_cellid_time_uid) cr->sort_by_cellid_time_uid();
			t1 = time(0);
			ctimec = ctime(&t1);
			for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
			fprintf(stderr,"%s, done sorting records\n",ctimec);
		}
	}
	
	// 2. load binary data if needed
	if(sfnin) {
		fprintf(stderr,"\treading serialized records from file %s\n",sfnin);
		
		// 2.1 open binary file, read data
		if(crecords::read_crecords_serialized(*cr2,sfnin,sinzip)) {
			fprintf(stderr,"Error reading records!\n");
			return 2;
		}
		t1 = time(0);
		ctimec = ctime(&t1);
		for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
		fprintf(stderr,"%s, %lu records read\n",ctimec,cr2->getnrecords());
				
		// 2.3 limit user IDs to the given ones if needed (this could be included in the records_read_tsv functionality as well)
		if(userlist || celllist || (tmin > 0 && tmax > tmin) ) {
			crecords* tmp = new crecords();
			tmp->copy_places_map(*cr2);
			
			filter_iterator<crecords::const_iterator,crecords::sentinel> fi(cr2->begin(),cr2->end());
			if(userlist) fi.set_userfilter(&userset);
			if(celllist) fi.set_cellsfilter(&cellset);
			if(tmin > 0 && tmax > tmin) fi.set_tsfilter(tmin,tmax);
			
			tmp->add_records_custom(fi,record_iterator_sentinel());
			
			delete cr2;
			cr2 = tmp;
		}
		
		// 2.4 sort data if needed
		if(sort_by_time || sort_by_uid_time || sort_by_cellid_time || sort_by_cellid_time_uid) {
			fprintf(stderr,"\tsorting records\n");
			if(sort_by_time) cr2->sort_by_time();
			if(sort_by_uid_time) cr2->sort_by_uid_time();
			if(sort_by_cellid_time) cr2->sort_by_cellid_time();
			if(sort_by_cellid_time_uid) cr2->sort_by_cellid_time_uid();
			t1 = time(0);
			ctimec = ctime(&t1);
			for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
			fprintf(stderr,"%s, done sorting serialized records\n",ctimec);
		}
	}
	
	
	// 3. create random data if needed
	if(randomtest) {
		if(sfnin) crrandom = cr;
		else crrandom = cr2;
		t1 = time(0);
		ctimec = ctime(&t1);
		for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
		fprintf(stderr,"%s, generating random data\n",ctimec);
		srand48(seed);
		crecords_fill_random_data(crrandom,rnum);
		
		t1 = time(0);
		ctimec = ctime(&t1);
		for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
		fprintf(stderr,"%s, generated %lu random records\n",ctimec,crrandom->getnrecords());
		
		// sort data if needed
		if(sort_by_time || sort_by_uid_time || sort_by_cellid_time || sort_by_cellid_time_uid) {
			fprintf(stderr,"\tsorting records\n");
			if(sort_by_time) crrandom->sort_by_time();
			if(sort_by_uid_time) crrandom->sort_by_uid_time();
			if(sort_by_cellid_time) crrandom->sort_by_cellid_time();
			if(sort_by_cellid_time_uid) crrandom->sort_by_cellid_time_uid();
			t1 = time(0);
			ctimec = ctime(&t1);
			for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
			fprintf(stderr,"%s, done sorting serialized records\n",ctimec);
		}
	}
	
	// 4. compare data (if more than one input was given)
	//	note: it does not make a lot of sense to compare with random data, but it can be done for testing purposes
	//	(if random output is saved and reloaded; it should be the save if run with the same seed both times)
	if( (ncdrfiles && sfnin) || (ncdrfiles && randomtest) || (sfnin && randomtest) ) {
		if(cr->compare_crecords(*cr2,cmp_ignore_places) == true) {
			fprintf(stdout,"OK\n");
		}
		else {
			if(sfnout) fprintf(stderr,"\nThere was an error comparing data; no data was saved\n");
			return 3; // make sure that if there was an error, output is not saved
		}
	}
	
	// 5. save binary data if required
	if(sfnout) {
		if(sfnin) crsave = cr2;
		else crsave = cr;
		
		//save crsave using the serialization
		if(crecords::write_crecords_serialized(*crsave,sfnout)) {
			fprintf(stderr,"Error writing records!\n");
			return 4;
		}
		t1 = time(0);
		ctimec = ctime(&t1);
		for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
		fprintf(stderr,"%s, done writing serialized records\n",ctimec);
	}
	
	delete cr; // note: these were always created (could be empty)
	delete cr2;
	
	return 0;
}



//compare function implemented here
//	note: if ignore_places == true, the places_map and all related fields are ignored
bool crecords::compare_crecords(crecords& other, bool ignore_places) {
			
	if(nrecords != other.nrecords) {
		fprintf(stderr,"crecords::compare_crecords(): nrecords differs (%lu != %lu)!\n",nrecords,other.nrecords);
		return false;
	}
	
	if(startstop.size() != other.startstop.size()) {
		fprintf(stderr,"crecords::compare_crecords(): size of start / end vectors differ (%lu != %lu)!\n",startstop.size(),other.startstop.size());
		return false;
	}
	
	if(nrecords != uids.size() || nrecords != ts.size() || nrecords != cells.size()) {
		fprintf(stderr,"crecords::compare_crecords(): nrecords differs from vector sizes (in this, %lu, %lu, %lu, %lu)!\n",
			nrecords,uids.size(),ts.size(),cells.size());
		return false;
	}
	
	if(other.nrecords != other.uids.size() || nrecords != other.ts.size() || nrecords != other.cells.size()) {
		fprintf(stderr,"crecords::compare_crecords(): nrecords differs from vector sizes (in other, %lu, %lu, %lu, %lu)!\n",
			other.nrecords,other.uids.size(),other.ts.size(),other.cells.size());
		return false;
	}
	
	for(uint64_t i=0;i<nrecords;i++) {
		if(uids[i] != other.uids[i]) {
			fprintf(stderr,"crecords::compare_crecords(): uids[%lu] differ: %ld != %ld!\n",i,uids[i],other.uids[i]);
			return false;
		}
		if(ts[i] != other.ts[i]) {
			fprintf(stderr,"crecords::compare_crecords(): ts[%lu] differ: %u != %u!\n",i,ts[i],other.ts[i]);
			return false;
		}
		if(cells[i] != other.cells[i]) {
			fprintf(stderr,"crecords::compare_crecords(): cells[%lu] differ: %d != %d!\n",i,cells[i],other.cells[i]);
			return false;
		}
	}
	
	if(sorted != other.sorted) {
		fprintf(stderr,"crecords::compare_crecords(): sort order differs!\n");
		return false;
	}
	if(uids_idx != other.uids_idx) {
		fprintf(stderr,"crecords::compare_crecords(): uids_idx differs!\n");
		return false;
	}
	if(cells_idx != other.cells_idx) {
		fprintf(stderr,"crecords::compare_crecords(): cells_idx differs!\n");
		return false;
	}
	if(uids_idx_cnts != other.uids_idx_cnts) {
		fprintf(stderr,"crecords::compare_crecords(): uids_idx_cnts differs!\n");
		return false;
	}
	
	/*
	 * variables to check for:
				uids,
				ts,
				cells,
				nrecords,
				
				sorted,
				uids_idx,
				cells_idx,
				uids_idx_cnts
	*/
	return true;
}

