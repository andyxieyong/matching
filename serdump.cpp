/*
 * serdump.cpp
 * 	convert back a serialized representation of records to plain text
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

int main(int argc, char **argv)
{
	char* sfn = 0; // input
	char* pfn = 0; // write out places_map separately if given
	bool inzip = false;
	
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'i':
			sfn = argv[i+1];
			break;
		case 'p':
			pfn = argv[i+1];
			break;
		case 'z':
			inzip = true;
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	
	crecords cr;
	if(crecords::read_crecords_serialized(cr,sfn,inzip)) {
		fprintf(stderr,"Error reading records!\n");
		return 2;
	}
	else fprintf(stderr,"%lu records read from file %s\n",cr.getnrecords(),sfn);
	
	if(pfn) {
		FILE* f = fopen(pfn,"w");
		if(!f) {
			fprintf(stderr,"Error opening output file %s!\n",pfn);
		}
		else {
			cr.print_places_map(f);
			fclose(f);
		}
	}
	
	record r;
	for(uint64_t i=0;cr.get_record(i,&r);i++) {
		if(r.startstop == 0 || r.startstop == 1) fprintf(stdout,"%ld\t%u\t%d\t%hhd\n",r.uid,r.ts,r.cell_id,r.startstop);
		else fprintf(stdout,"%ld\t%u\t%d\n",r.uid,r.ts,r.cell_id);
	}
	
	return 0;
}

