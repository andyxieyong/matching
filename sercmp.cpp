/*
 * sercmp.cpp -- compare two binary (serialized) archives
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
#include <cereal/archives/portable_binary.hpp>


int main(int argc, char **argv)
{
	char* sfn1 = 0; //file name for serialized data (first)
	char* sfn2 = 0; //file name for the second
	bool zip1 = false; //first file is gzipped
	bool zip2 = false; //second file is gzipped
	
	// perform explicit sort on both inputs before comparing
	bool sort_by_uid_time = false; //sort the input by uid and time
	bool sort_by_time = false; //sort the input by time
	bool sort_by_cellid_time = false; //sort the input by cellid and time
	bool sort_by_cellid_time_uid = false; //sort by cellid and time and also uid (so that the sorted result should be unique)
	
	char* ctimec = 0;
	time_t t1;
	
	struct popen_noshell_pass_to_pclose spc;
	char* popen_args[4];
	popen_args[0] = "/bin/gzip";
	popen_args[1] = "-cd";
	popen_args[2] = 0;
	popen_args[3] = 0;
	
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'z':
		case 'Z':
			if(argv[i][2] == '1') { zip1 = true; break; }
			if(argv[i][2] == '2') { zip2 = true; break; }
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
		case '1':
			sfn1 = argv[i+1];
			i++;
			break;
		case '2':
			sfn2 = argv[i+1];
			i++;
			break;
		case 'S':
			if(argv[i][2] == 'u') { sort_by_uid_time = true; break; }
			if(argv[i][2] == 't') { sort_by_time = true; break; }
			if(argv[i][2] == 'c') { sort_by_cellid_time = true; break; }
			if(argv[i][2] == 'C') { sort_by_cellid_time_uid = true; break; }
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	
	if(sfn1 == 0 || sfn2 == 0) {
		fprintf(stderr,"No binary file name specified!\n");
		return 1;
	}
	
	// 1. load from the two binary files
	crecords cr;
	crecords cr2;
	
	//load previously saved instance, compare to cr
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"\%s, reading serialized records from file %s\n",ctimec,sfn1);
	
	if(zip1) {
		popen_args[2] = sfn1;
		FILE* f = popen_noshell(popen_args[0],popen_args,"r",&spc,0);
		stdiostream sifs(f);
		cereal::PortableBinaryInputArchive ar(sifs.stream());
		ar(cr);
		pclose_noshell(&spc); popen_args[2] = 0;
	}
	else {
		std::ifstream ifs(sfn1,std::ifstream::binary);
		cereal::PortableBinaryInputArchive ar(ifs);
		ar(cr);
	}
	
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"\%s, %lu records read\n\treading serialized records from file %s\n",ctimec,cr.getnrecords(),sfn2);
	
	if(zip2) {
		popen_args[2] = sfn2;
		FILE* f = popen_noshell(popen_args[0],popen_args,"r",&spc,0);
		stdiostream sifs(f);
		cereal::PortableBinaryInputArchive ar(sifs.stream());
		ar(cr2);
		pclose_noshell(&spc); popen_args[2] = 0;
	}
	else {
		std::ifstream ifs(sfn2,std::ifstream::binary);
		cereal::PortableBinaryInputArchive ar(ifs);
		ar(cr2);
	}
	
	t1 = time(0);
	ctimec = ctime(&t1);
	for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
	fprintf(stderr,"%s, %lu records read\n",ctimec,cr2.getnrecords());
	
	// sort the input if needed before comparing
	if(sort_by_time || sort_by_uid_time || sort_by_cellid_time || sort_by_cellid_time_uid) {
		fprintf(stderr,"\tsorting records\n");
		cr.use_combined_sort = true;
		if(sort_by_time) cr.sort_by_time();
		if(sort_by_uid_time) cr.sort_by_uid_time();
		if(sort_by_cellid_time) cr.sort_by_cellid_time();
		if(sort_by_cellid_time_uid) cr.sort_by_cellid_time_uid();
		
		cr2.use_combined_sort = true;
		if(sort_by_time) cr2.sort_by_time();
		if(sort_by_uid_time) cr2.sort_by_uid_time();
		if(sort_by_cellid_time) cr2.sort_by_cellid_time();
		if(sort_by_cellid_time_uid) cr2.sort_by_cellid_time_uid();
		
		t1 = time(0);
		ctimec = ctime(&t1);
		for(char* c1 = ctimec;*c1;c1++) if(*c1 == '\n') *c1 = 0;
		fprintf(stderr,"%s, done sorting records\n",ctimec);
	}
	
	if(cr.compare_crecords(cr2) == true) {
		fprintf(stdout,"OK\n");
	}
	
	
	return 0;
}



//compare function implemented here
bool crecords::compare_crecords(crecords& other) {
	
	if(nrecords != other.nrecords) {
		fprintf(stderr,"crecords::compare_crecords(): nrecords differs (%lu != %lu)!\n",nrecords,other.nrecords);
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
	
	if(places_map != other.places_map) {
		fprintf(stderr,"crecords::compare_crecords(): places_map differs!\n");
		return false;
	}
	
	if(cell_coords != other.cell_coords) {
		fprintf(stderr,"crecords::compare_crecords(): cell_coords differs!\n");
		return false;
	}
	
	if(vlon != other.vlon) {
		fprintf(stderr,"crecords::compare_crecords(): vlon differs!\n");
		return false;
	}
	
	if(vlat != other.vlat) {
		fprintf(stderr,"crecords::compare_crecords(): vlat differs!\n");
		return false;
	}
	if(vertex_cell != other.vertex_cell) {
		fprintf(stderr,"crecords::compare_crecords(): vertex_cell differs!\n");
		return false;
	}
	if(cell_dup != other.cell_dup) {
		fprintf(stderr,"crecords::compare_crecords(): cell_dup differs!\n");
		return false;
	}
	if(cell_dup_rev != other.cell_dup_rev) {
		fprintf(stderr,"crecords::compare_crecords(): cell_dup_rev differs!\n");
		return false;
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
				
				places_map,
				cell_coords,
				vlon,
				vlat,
				vertex_cell,
				edges_cell,
				cell_dup,
				cell_dup_rev,
				
				sorted,
				uids_idx,
				cells_idx,
				uids_idx_cnts
	*/
	return true;
}

