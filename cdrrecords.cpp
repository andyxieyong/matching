/*
 * cdrrecords.cpp
 * 
 * read CDR or similar checkin records, perform spatiotemporal search
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
#include <string.h>
#include <ctype.h>
#include <stdlib.h>



// search in time (binary search, limit between s and e)
uint64_t crecords::tsearch(uint32_t t1, bool start, uint64_t s, uint64_t e) const {
	// use std implementations for search instead custom
	
	// binary search in the ts array (note: array needs to be sorted by time, this is currently ensured by the read functions
	//		and the functions calling this one)
	// return value:
	//		if start == true, the first occurrence of t1 or the smallest ts value > t1 if it is not found
	//			or e if t1 > all ts values
	//		if start == false, the last occurrence of t1 or the largest ts value < t1 if it is not found
	//			or e if t1 < all ts values (note: 0 would mean that ts[0] <= t1)
	// this was a range can be found with searching for tsearch(start,true) and tsearch(end,false);
	// limit search between s and e (s inclusive, e exclusive; if omitted, search the whole range)
	
	// std::lower_bound:  Returns an iterator pointing to the first element in the range [first,last) which does not compare less than val.
	//	in the case of start == true, this is the return value
	// std::upper_bound: Returns an iterator pointing to the first element in the range [first,last) which compares greater than val.
	//	in the case of start == false, this is the value coming after the return value
	
	if(s >= ts.size()) return e; // no records to search in (note: this also covers when ts.size() == 0, since s >= 0 always)
	if(e > ts.size()) e = ts.size();
	
	if(start) {
		std::vector<uint32_t>::const_iterator it = std::lower_bound(ts.begin() + s, ts.begin() + e, t1);
		return it - ts.begin();
	}
	else {
		if(ts[s] > t1) return e; // otherwise, we would end up returning 0-1 in the next lines
		std::vector<uint32_t>::const_iterator it = std::upper_bound(ts.begin() + s, ts.begin() + e, t1);
		return (it - ts.begin()) - 1; // note: the caller expects the last element <= t1, so we subtract 1; the case when it == ts.start() is taken care of above
	}
}


// copy a places_map, reversing mappings
static void reverse_places_map(hashmap<int,hashset<int> >& new_map, const hashmap<int,hashset<int> >& old_map) {
	for(hashmap<int,hashset<int> >::const_iterator it = old_map.begin(); it != old_map.end(); it++) {
		int pid = it->first;
		const hashset<int>& cells = it->second;
		for(hashset<int>::const_iterator it2 = cells.begin(); it2 != cells.end(); it2++)
			new_map[*it2].insert(pid); //note: duplicates are taken care of automatically, and a new, empty hashset is also created if needed
	}
}


// read a place-cell mapping from the specified file and add it to the list (places_maps)
// return values: 0: success; >0: error
// it is considered an error if the type is different than what was previously read or if the radius is exactly the same as previously
int crecords::read_places_map(FILE* in) {
	// first two lines: header: type and radius
	char* line1 = 0;
	char* line2 = 0;
	size_t n = 0;
	
	if(getline(&line1,&n,in) == -1) return 1;
	if(line1 == 0) return 1;
	line2 = line1;
	while(isspace(*line2)) line2++;
	if(*line2 == 0) { free(line1); return 1; } // empty first line, this is an error
	if(*line2 != '#') { free(line1); return 1; } // format: # (comment sign), type description
	do line2++; while(isspace(*line2)); // skip '#' and further spaces
	
	// determine type
	places_map_type pmt1;
	if( ! strncmp(line2,"place cell",10) ) { pmt1 = place_cell; goto cmp_end; }
	if( ! strncmp(line2,"place place",11) ) { pmt1 = place_place; goto cmp_end; }
	if( ! strncmp(line2,"cell cell",9) ) { pmt1 = cell_cell; goto cmp_end; }
	if( ! strncmp(line2,"cell cell (centers)",19) ) { pmt1 = cell_cell_centers; goto cmp_end; }
	return 2; // if we are here, the format is invalid
	
cmp_end:
	if( places_maps.size() > 0 && pmt1 != pmt ) {
		fprintf(stderr,"crecords::read_places_map(): incompatible map type in the given file!\n");
		free(line1);
		return 3;
	}
	else pmt = pmt1;
	
	if(getline(&line1,&n,in) == -1) { free(line1); return 1; }
	if(line1 == 0) return 1;
	line2 = line1;
	while(isspace(*line2)) line2++;
	if(*line2 == 0) { free(line1); return 1; } // empty first line, this is an error
	if(*line2 != '#') { free(line1); return 1; } // format: # (comment sign), type description
	do line2++; while(isspace(*line2)); // skip '#' and further spaces
	
	// read radius
	double radius;
	if(sscanf(line2,"%lf",&radius) != 1) { free(line1); return 2; }
	
	// check if this radius is not alread stored
	for(unsigned int i=0;i<radii.size();i++) if(radii[i] == radius) {
		fprintf(stderr,"crecords::read_places_map(): radius (%g) already present!\n",radius);
		free(line1);
		return 3;
	}
	
	radii.push_back(radius);
	places_maps.push_back(hashmap<int,hashset<int> >());
	hashmap<int,hashset<int> >& places_map = places_maps.back();
		
	unsigned int l = 2;
	while(1) {
		if(getline(&line1,&n,in) == -1) break;
		if(line1 == 0) break;
		l++;
		int c1,c2;
		if(sscanf(line1,"%d %d",&c1,&c2) != 2) {
			fprintf(stderr,"crecords::read_places_map(): invalid data in line %u!\n",l);
			free(line1);
			return 4;
		}
		places_map[c1].insert(c2);
	}
	if(line1) free(line1);
	fprintf(stderr,"\tcrecords:read_places_map(): %u records read\n",l-2);
	if(places_map_default == 0) places_map_default = &(places_maps[0]); // make sure the default is set
	return 0;
}
		

// print out matching between places and cells -- note: reading back will not work for files containing more than one map
//	also for the cell -- place type
void crecords::print_places_map(FILE* out) const {
	for(unsigned int i=0;i<radii.size();i++) {
		switch(pmt) {
			case place_cell: fprintf(out,"# place cell\n"); break;
			case place_place: fprintf(out,"# place place\n"); break;
			case cell_cell: fprintf(out,"# cell cell\n"); break;
			case cell_cell_centers: fprintf(out,"# cell cell (centers)\n"); break;
			case cell_place: fprintf(out,"# cell place\n"); break;
		}
		fprintf(out,"# %g\n",radii[i]);
		const hashmap<int,hashset<int> >& places_map = places_maps[i];
		for(hashmap<int,hashset<int> >::const_iterator it = places_map.begin(); it != places_map.end(); it++) {
			int pid = it->first;
			const hashset<int>& cells = it->second;
			for(hashset<int>::const_iterator it2 = cells.begin(); it2 != cells.end(); it2++)
				fprintf(out,"%d\t%d\n",pid,*it2);
		}
	}
}


void crecords::copy_places_map_reverse(const crecords& cr) {
	places_maps.clear();
	radii.clear();
	
	for(unsigned int i=0;i<cr.radii.size();i++) {
		radii.push_back(cr.radii[i]);
		places_maps.push_back(hashmap<int,hashset<int> >());
		hashmap<int,hashset<int> >& places_map = places_maps.back();
		
		reverse_places_map(places_map, cr.places_maps[i]);
	}
	
	if(cr.pmt == place_cell || cr.pmt == cell_place) {
		if(cr.pmt == place_cell) pmt = cell_place;
		else pmt = place_cell;
	}
	else pmt = cr.pmt;
}


 //fill up cell_idx, so that searching by cell_id is made easier
void crecords::create_cell_idx(){
	if( ! (sorted == cellid_time || sorted == cellid_time_uid) ) {
		fprintf(stderr,"crecords::create_cell_idx(): records are not properly sorted!\n");
		return;
	}
	if(cells.size() == 0) return;
	uint64_t i = 0;
	uint64_t j = 1;
	for(;;j++) {
		if(j == cells.size()) {
			cells_idx[ cells[i] ] = std::make_pair( i, j-i );
			break;
		}
		if(cells[j] != cells[i]) {
			cells_idx[ cells[i] ] = std::make_pair( i, j-i );
			i = j;
		}
	}
}


 //fill up uids_idx, so that searching by uid is made easier
void crecords::create_uid_idx(){
	if( ! ( sorted == uid || sorted == uid_time || sorted == uid_cellid || sorted == uid_time_cellid ) ) {
		fprintf(stderr,"crecords::create_uid_idx(): records are not properly sorted!\n");
		return;
	}
	uids_idx.clear();
	uids_idx_cnts = false;
	if(uids.size() == 0) return;
	uint64_t i = 0;
	uint64_t j = 1;
	for(;;j++) {
		if(j == uids.size()) {
			uids_idx[ uids[i] ] = std::make_pair( i, j-i );
			break;
		}
		if(uids[j] != uids[i]) {
			uids_idx[ uids[i] ] = std::make_pair( i, j-i );
			i = j;
		}
	}
}



// set the parameters used for comparisons by the above functions
// need to be called before calling any of those, but after the coordinates and the Voronoi-tesselation has been loaded
// not thread-safe (since it modifies the stored params), should be called by the main thread at first
// generate appropriate places_maps if needed
// if searchreverse == true, the mappings are created reversed (cell_id -> place_id), this should be used if this class contains transportation records
int crecords::set_compare_params(unsigned int dt1_, unsigned int dt2_, double radius1p, double radius1i, double radius2p, double radius2i, bool searchreverse) {
	places_map_1p = 0; places_map_1i = 0; places_map_2p = 0; places_map_2i = 0;
	if( ! (radius1p == 0.0 && radius1i == 0.0 && radius2p == 0.0 && radius2i == 0.0) ) {
	
		if(places_maps.size() == 0) {
			fprintf(stderr,"crecords::set_compare_params(): error: no place-cell mappings loaded!\n");
			return -1;
		}
		if( (searchreverse == false && pmt != place_cell) || (searchreverse == true && pmt != cell_place) ) {
			fprintf(stderr,"crecords::set_compare_params(): error: incompatible map types!\n");
			return -1;
		}
		
		
		// find the approrpiate maps and store them
		for(unsigned int i=0;i<radii.size();i++) {
			if(radii[i] == radius1p) places_map_1p = &(places_maps[i]);
			if(radii[i] == radius1i) places_map_1i = &(places_maps[i]);
			if(radii[i] == radius2p) places_map_2p = &(places_maps[i]);
			if(radii[i] == radius2i) places_map_2i = &(places_maps[i]);
			
			if( places_map_1p && places_map_1i && places_map_2p && places_map_2i ) break; // all has been assigned already
		}
	
		if( places_map_1p == 0 ) fprintf(stderr,"crecords::set_compare_params(): error: radius %g missing!\n",radius1p);
		if( places_map_1i == 0 ) fprintf(stderr,"crecords::set_compare_params(): error: radius %g missing!\n",radius1i);
		if( places_map_2p == 0 ) fprintf(stderr,"crecords::set_compare_params(): error: radius %g missing!\n",radius2p);
		if( places_map_2i == 0 ) fprintf(stderr,"crecords::set_compare_params(): error: radius %g missing!\n",radius2i);
		if( !( places_map_1p && places_map_1i && places_map_2p && places_map_2i ) ) return -1;
	}
	
	if(radius1p > radius2p) places_map_default = places_map_1p;
	else places_map_default = places_map_2p;
	
	dt1 = dt1_; // time period for "walking": used before a record if startstop == false, used after a record when startstop == true
	dt2 = dt2_; // time period for "transit": used before a record if startstop == true, used after a record when startstop == false
	return 0;	
}

// set only the default map
int crecords::set_default_places_map(double radius) {
	if(places_maps.size() == 0) {
		fprintf(stderr,"crecords::set_default_places_map(): error: no place-cell mappings loaded!\n");
		return -1;
	}
	
	places_map_default = 0;
	for(unsigned int i=0;i<radii.size();i++) {
		if(radii[i] == radius) {
			places_map_default = &(places_maps[i]);
			break;
		}
	}
	
	if(places_map_default == 0) {
		fprintf(stderr,"crecords::set_default_places_map(): error: radius %g missing!\n",radius);
		return -1;
	}
	return 0;
}


// count all matches for events
// number of possible matches returned in *matches_possible
// if only_temporal is true, spatial consistency is not checked and *matches_possible contains all temporal matches
// if count_impossible == false (default), than one impossible match aborts the compaision and 0 is returned (i.e. matches_possible and
//		matches_impossible are not updated); if all matches are possible, then *matches_possible contains the number of possible (temporal)
//		matches, but matches_impossible is not dereferenced (i.e. is count_impossible == false, then matches_impossible can be NULL or invalid)
// if count_impossible == true, and there is at least one impossible match then *matches_possible contains the number of temporal matches (as if
//		only_temporal == true was given); if matches_impossible is not NULL, than the number of impossible matches is counted and stored there
// return value: 1: possible (*matches is the number of matches), 0: impossible (at least one inconsistent match), -1: error while searching (place_id not found)
int crecords::compare_users2(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches_possible,
		unsigned int* matches_impossible, bool only_temporal, bool ignore_missing, bool count_impossible, uint32_t tmin, uint32_t tmax) const {
	
	if(tmin > tmax) {
		fprintf(stderr,"crecords::compare_users2(): error: tmin > tmax!\n");
		return -1;
	}
	
	//note: both records must be sorted by uid and time
	if( ! (sorted == uid_time || sorted == uid_time_cellid) ) {
		fprintf(stderr,"crecords::compare_users2(): error: records not properly sorted!\n");
		return -1;
	}
	if( ! (other_records.sorted == uid_time || other_records.sorted == uid_time_cellid) ) {
		fprintf(stderr,"crecords::compare_users2(): error: other records not properly sorted!\n");
		return -1;
	}
	
	//find uids in both records
	hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it1 = uids_idx.find(uid);
	if(it1 == uids_idx.end()) {
		fprintf(stderr,"crecords::compare_users2(): error: uid %ld not found!\n",uid);
		return -1;
	}
	hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it2 = other_records.uids_idx.find(other_uid);
	if(it2 == other_records.uids_idx.end()) {
		fprintf(stderr,"crecords::compare_users2(): error: other uid %ld not found!\n",other_uid);
		return -1;
	}
	
	if( !( (places_map_1i && places_map_1p && places_map_2i && places_map_2p) || only_temporal) ) {
		fprintf(stderr,"crecords::compare_users2(): places_maps not set up properly (use the set_compare_params() function before calling this)!\n");
		return -1;
	}
	
	// go through all records (in a possibly not so efficient manner, i.e. do a binary search for each event instead of iterating through the two
	//	vectors in parallel
	uint64_t s = it2->second.first;
	uint64_t e = it2->second.first + it2->second.second;
	uint64_t s2 = it1->second.first;
	uint64_t e2 = it1->second.first + it1->second.second;
	
	
	// adjust first and last record if the time interval is limited
	if(tmin > 0 && tmax > 0) {
		s = other_records.tsearch(tmin,true,s,e);
		if(s < e) e = other_records.tsearch(tmax,false,s,e) + 1; // + 1, since tsearch returns the result inclusive, but expects the end of the range exclusive
		// adjustment: ts[e-1] == tmax is possible, this needs to be excluded
		while(e > s && other_records.ts[e-1] == tmax) e--;
		if(s >= e) { // nothing to compare, the other user has no records in this time frame
			*matches_possible = 0;
			if(count_impossible && matches_impossible) *matches_impossible = 0;
			return 1;
		}
	}
	
	if(tmin > 0 && tmax > 0) {
		s2 = tsearch(tmin,true,s2,e2);
		if(s2 < e2) e2 = tsearch(tmax,false,s2,e2) + 1;
		// adjustment: ts[e-1] == tmax is possible, this needs to be excluded
		while(e2 > s2 && ts[e2-1] == tmax) e2--;
		if(s2 >= e2) {
			*matches_possible = 0;
			if(count_impossible && matches_impossible) *matches_impossible = 0;
			return 1;
		}
	}
	
	unsigned int mp = 0;
	unsigned int mi = 0;
	int r = 1;
	for(uint64_t i = s2; i < e2; i++) {
		unsigned int ts11 = ts[i] - dt1;
		unsigned int ts12 = ts[i] + dt1;
		
		// if trip start / stop info is used, apply the separate time period for the "transit" trips
		if(startstop.size() > 0) {
			if(startstop[i]) ts11 = ts[i] - dt2; // end of the trip, transit time before the current ts
			else ts12 = ts[i] + dt2; // start of the trip, transit time after the current ts
		}
		
		// also limit the time period according to what is given as parameters
		if(tmin > 0 && tmax > 0) {
			if(ts11 < tmin) ts11 = tmin;
			if(ts12 > tmax) ts12 = tmax;
		}
		
		// check the previous and next records and limit the time period to not go after those
		if(i > s2) if(ts11 <= ts[i-1]) ts11 = ts[i-1] + 1;
		if(i < e2-1) if(ts12 >= ts[i+1]) ts12 = ts[i+1] - 1;
		
		uint64_t start = other_records.tsearch(ts11, true, s, e);
		if(start >= e) continue; //not found (no CDR events for this user in this time window)
		uint64_t end = other_records.tsearch(ts12, false, s, e);
		if(end >= e) continue;
		
		// at this point, there is at least one temporal match
		mp++;
		if(only_temporal) continue; // in this case, we do not need to check spatial compatibility
		
		bool possible_found = false;
		bool impossible_found = false;
		
		const hashmap<int, hashset<int> >* map_possible = places_map_1p;
		const hashmap<int, hashset<int> >* map_impossible = places_map_1i;
		
		
		if(startstop.size() > 0) if(startstop[i]) {
			map_possible = places_map_2p;
			map_impossible = places_map_2i;
		}
		
		for(uint64_t j = start; j <= end; j++) {
			if(other_records.ts[j] > ts[i] && startstop.size() > 0) {
				if(startstop[i]) {
					map_possible = places_map_1p;
					map_impossible = places_map_1i;
				}
				else {
					map_possible = places_map_2p;
					map_impossible = places_map_2i;
				}
			}
			
			hashmap<int, hashset<int> >::const_iterator it3 = map_possible->find(other_records.cells[j]);
			if(it3 != map_possible->end())
				if(it3->second.count(cells[i]) > 0) { possible_found = true; continue; }
			
			hashmap<int, hashset<int> >::const_iterator it4 = map_impossible->find(other_records.cells[j]);
			if(it4 == map_impossible->end()) {
				if(ignore_missing) { impossible_found = true; break; } //consider this as impossible here too
				fprintf(stderr,"crecords::compare_users(): error: place_id %d not found!\n",other_records.cells[j]);
				return -1;
			}
				
			if(it4->second.count( cells[i] ) == 0)
				{ impossible_found = true; break; } //impossible (note: in this case, there is no need to check the others as well)
		}
		if(impossible_found) {
			if(count_impossible) {
				if(matches_impossible) mi++;
				else only_temporal = true;
				r = 0;
			}
			else return 0;
		}
		
	}
	
	*matches_possible = mp;
	if(count_impossible && matches_impossible) *matches_impossible = mi;
	return r; // in this case, we only have possible matches
}



// compare record i in this with record other_i in other_records, using the parameters set in this
// return value:
//	1: other_i is a match
//	0: other_i is not a temporal match or neither possible nor impossible (i.e. not in places_map_?p, but in places_map_?i)
//	-1: other_i is an impossible match (temporal match but not spatial)
//	-2: error (place_id not found)
// if only_temporal == true: 1 if it is a temporal match, 0 otherwise
int crecords::compare_records(uint64_t i, const crecords& other_records, uint64_t other_i, bool ignore_missing, bool only_temporal) const {
	const hashmap<int, hashset<int> >* map_possible = places_map_1p;
	const hashmap<int, hashset<int> >* map_impossible = places_map_1i;
	unsigned int dt = dt1;
	unsigned int other_ts = other_records.ts[ other_i ];
	
	if(startstop.size() > 0) {
		if( 
				( startstop[i] == true && ts[i] >= other_ts ) // end of the trip, other record was during the trip
			||	( startstop[i] == false && ts[i] < other_ts ) // start of the trip, other record is after
			) { // transit instead of walking
				dt = dt2;
				map_possible = places_map_2p;
				map_impossible = places_map_2i;
			}
	}
	
	// check if it is really a temporal match (note: dt already corresponds to whether ts[i] is greater or smaller than other_ts)
	if( ts[i] + dt < other_ts || ts[i] - dt > other_ts ) return 0;
	
	// check if there is an overlap with the next / previous record
	if( timeoverlap == false ) {
		if( ts[i] < other_ts ) if(i+1 < nrecords && uids[i+1] == uids[i]) if(ts[i+1] <= other_ts) return 0;
		if( ts[i] > other_ts ) if(i > 0 && uids[i-1] == uids[i]) if(ts[i-1] >= other_ts) return 0;
	}
	if( other_records.timeoverlap == false ) {
		if( ts[i] < other_ts ) if(other_i > 0 && other_records.uids[other_i-1] == other_records.uids[other_i])
			if(other_records.ts[other_i-1] >= ts[i]) return 0;
		if( ts[i] > other_ts ) if(other_i+1 < other_records.nrecords && other_records.uids[other_i+1] == other_records.uids[other_i])
			if(other_records.ts[other_i+1] <= ts[i]) return 0;
	}
	
	if(only_temporal) return 1; // this is now a certain temporal match
	
	// check spatial match
	int other_cellid = other_records.cells[other_i];
	hashmap<int, hashset<int> >::const_iterator it = map_possible->find( other_cellid );
	if(it != map_possible->end()) if(it->second.count( cells[i] ) > 0) return 1; // possible match
	
	hashmap<int, hashset<int> >::const_iterator it2 = map_impossible->find( other_cellid );
	if(it2 == map_impossible->end()) {
		if(ignore_missing) return -1; // treat as impossible as well
		fprintf(stderr,"crecords::compare_records(): error: place_id %d not found!\n",other_cellid);
		return -2;
	}
	
	if(it2->second.count(cells[i]) > 0) return 0; //undecided
	else return -1; // impossible
}


// count matches, but count each event from each dataset only once maximum
// number of possible matches returned in *matches_possible
// if only_temporal is true, spatial consistency is not checked and *matches_possible contains all temporal matches
// if count_impossible == false (default), than one impossible match aborts the compaision and 0 is returned (i.e. matches_possible and
//		matches_impossible are not updated); if all matches are possible, then *matches_possible contains the number of possible (temporal)
//		matches, but matches_impossible is not dereferenced (i.e. is count_impossible == false, then matches_impossible can be NULL or invalid)
// if count_impossible == true, and there is at least one impossible match then *matches_possible contains the number of temporal matches (as if
//		only_temporal == true was given); if matches_impossible is not NULL, than the number of impossible matches is counted and stored there
// return value: 1: possible (*matches is the number of matches), 0: impossible (at least one inconsistent match), -1: error while searching (place_id not found)
int crecords::compare_users2_onematch(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches_possible,
		unsigned int* matches_impossible, bool only_temporal, bool ignore_missing, bool count_impossible, uint32_t tmin, uint32_t tmax) const {
	
	
	//note: both records must be sorted by uid and time
	if( ! (sorted == uid_time || sorted == uid_time_cellid) ) {
		fprintf(stderr,"crecords::compare_users_onematch(): error: records not properly sorted!\n");
		return -1;
	}
	
	if( ! (other_records.sorted == uid_time || other_records.sorted == uid_time_cellid) ) {
		fprintf(stderr,"crecords::compare_users_onematch(): error: other records not properly sorted!\n");
		return -1;
	}
	
	//find uids in both records
	hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it1 = uids_idx.find(uid);
	if(it1 == uids_idx.end()) {
		fprintf(stderr,"crecords::compare_users_onematch(): error: uid %ld not found!\n",uid);
		return -1;
	}
	hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it2 = other_records.uids_idx.find(other_uid);
	if(it2 == other_records.uids_idx.end()) {
		fprintf(stderr,"crecords::compare_users_onematch(): error: other uid %ld not found!\n",other_uid);
		return -1;
	}
	
	if( ! ((places_map_1i && places_map_1p && places_map_2i && places_map_2p) || only_temporal) ) {
		fprintf(stderr,"crecords::compare_users_onematch(): places_maps not set up properly (use the set_compare_params() function before calling this)!\n");
		return -1;
	}
	
	unsigned int mp = 0;
	unsigned int mi = 0;
	
	uint64_t s2 = it2->second.first;
	uint64_t e2 = it2->second.first + it2->second.second;
	uint64_t s1 = it1->second.first;
	uint64_t e1 = it1->second.first + it1->second.second;
	
	if(tmin > 0 && tmax > 0) {
		s2 = other_records.tsearch(tmin,true,s2,e2);
		if(s2 < e2) e2 = other_records.tsearch(tmax,false,s2,e2) + 1; // + 1, since tsearch returns the result inclusive, but expects the end of the range exclusive
		// adjustment: ts[e-1] == tmax is possible, this needs to be excluded
		while(e2 > s2 && other_records.ts[e2-1] == tmax) e2--;
		if(s2 >= e2) { // nothing to compare, the other user has no records in this time frame
			*matches_possible = 0;
			if(count_impossible && matches_impossible) *matches_impossible = 0;
			return 1;
		}
		s1 = tsearch(tmin,true,s1,e1);
		if(s1 < e1) e1 = tsearch(tmax,false,s1,e1) + 1;
		// adjustment: ts[e-1] == tmax is possible, this needs to be excluded
		while(e1 > s1 && ts[e1-1] == tmax) e1--;
		if(s1 >= e1) {
			*matches_possible = 0;
			if(count_impossible && matches_impossible) *matches_impossible = 0;
			return 1;
		}
	}
	
	unsigned int dt = dt1;
	if(dt2 > dt) dt = dt2;
	int r = 1; // return value (possible by default, adjusted if impossible matches are found)
	
	while(s1 < e1 && s2 < e2) {
		//current records: s1 in this and s2 in other_records
		
		bool found_possible = false;
		bool found_impossible = false;
		int r1;
		
		//do the processing for the first one (with respect to temporal order)
		if(ts[s1] < other_records.ts[s2]) {
			if(ts[s1] + dt < other_records.ts[s2]) {
				s1++;
				continue;
			}
			
			r1 = compare_records(s1,other_records,s2,ignore_missing,only_temporal);
			if(r1 == -2) return -1; // error
			if(only_temporal) {
				if(r1) { mp++; s2++; }
				s1++;
			}
			else {
				// we have to compare s1 to each record after s2 as well to check for impossible matches
				if(r1 == -1) found_impossible = true;
				else for(uint64_t s22 = s2+1; s22 < e2 && other_records.ts[s22] <= ts[s1] + dt; s22++) {
					int r2 = compare_records(s1,other_records,s22,ignore_missing,false);
					if(r2 == -2) return -1; // error
					if(r2 == -1) { found_impossible = true; break; }
				}
				
				if(r1 == 1 && !found_impossible) {
					found_possible = true;
					for(uint64_t s12 = s1+1; s12 < e1 && ts[s12] <= other_records.ts[s2] + dt; s12++) {
						int r2 = compare_records(s12,other_records,s2,ignore_missing,false);
						if(r2 == -2) return -1; // error
						if(r2 == -1) { found_impossible = true; break; }
					}
				}
				if(r1 == -1 || r1 == 1) s2++; // if s1 -- s2 is a temporal match, increase s2 as well
				s1++;
			}
		} //!! 
		else {
			if(ts[s1] > other_records.ts[s2] + dt) {
				s2++;
				continue;
			}
			
			r1 = compare_records(s1,other_records,s2,ignore_missing,only_temporal);
			if(r1 == -2) return -1; // error
			if(only_temporal) {
				if(r1) { mp++; s1++; }
				s2++;
			}
			else {
				if(r1 == -1) found_impossible = true;
				else for(uint64_t s12 = s1+1; s12 < e1 && ts[s12] <= other_records.ts[s2] + dt; s12++) {
					int r2 = compare_records(s12,other_records,s2,ignore_missing,false);
					if(r2 == -2) return -1; // error
					if(r2 == -1) { found_impossible = true; break; }
				}
				
				if(r1 == 1 && !found_impossible) {
					found_possible = true;
					for(uint64_t s22 = s2+1; s22 < e2 && other_records.ts[s22] <= ts[s1] + dt; s22++) {
						int r2 = compare_records(s1,other_records,s22,ignore_missing,false);
						if(r2 == -2) return -1; // error
						if(r2 == -1) { found_impossible = true; break; }
					}
				}
				if(r1 == -1 || r1 == 1) s1++; // if s1 -- s2 is a temporal match, increase s1 as well
				s2++;
				
				
			}
		}
		
		if(found_impossible) {
			if(count_impossible) {
				if(r1 != 0) mp++; // s1 -- s2 was a temporal match -> mp always counts the number of temporal matches in this case
				if(matches_impossible) mi++; // in this case, we exactly count the number of impossible matches
				else only_temporal = true; // otherwise, just proceed with counting temporal matches
				r = 0; // record that there was an impossible match
			}
			else return 0;
		}
		else if(found_possible) mp++;
	}
	
	
	*matches_possible = mp;
	if(count_impossible && matches_impossible) *matches_impossible = mi;
	return r;
}





// count matches, but count each event from each dataset only once maximum
// similarly to the previous version but also consider a (growing) sequence of time intervals,
//	given by the tmax array (length of ntmax)
// matches_temporal[i] stores the number of temporal matches between [tmin, tmax[i]), matches_possible[i] stores the
//	number of possible matches in this interval (or 0, if it is impossible)
// if only_temporal is true, spatial consistency is not checked and elements in matches_possible are not accessed at all
// return value: first tmax that there is at least one impossible match (i.e. if it is r, than there is at least one impossible match between
//	tmin and tmax[r]; if there are no impossible matches, it is ntmax), or -1 if there is an error processing
// note: it is an error if values in tmax are not ordered by time or are closer to each other than the larger of dt1 and dt2
// tmin and the last value in tmax can be 0, meaning the whole interval
int crecords::compare_users_onematch_tmulti(int64_t uid, const crecords& other_records, int64_t other_uid, int ntmax,
			unsigned int* matches_temporal, unsigned int* matches_possible, uint32_t tmin, const uint32_t* tmax,
			bool only_temporal, bool ignore_missing) const {
	
	//note: both records must be sorted by uid and time
	if( ! (sorted == uid_time || sorted == uid_time_cellid) ) {
		fprintf(stderr,"crecords::compare_users_onematch_tmulti(): error: records not properly sorted!\n");
		return -1;
	}
	
	if( ! (other_records.sorted == uid_time || other_records.sorted == uid_time_cellid) ) {
		fprintf(stderr,"crecords::compare_users_onematch_tmulti(): error: other records not properly sorted!\n");
		return -1;
	}
	
	if(ntmax <= 0 || tmax == 0 || matches_temporal == 0 || (only_temporal == false && matches_possible == 0) ) {
		fprintf(stderr,"crecords::compare_users_onematch_tmulti(): invalid parameters!\n");
		return -1;
	}
	
	//find uids in both records
	hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it1 = uids_idx.find(uid);
	if(it1 == uids_idx.end()) {
		fprintf(stderr,"crecords::compare_users_onematch_tmulti(): error: uid %ld not found!\n",uid);
		return -1;
	}
	hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it2 = other_records.uids_idx.find(other_uid);
	if(it2 == other_records.uids_idx.end()) {
		fprintf(stderr,"crecords::compare_users_onematch_tmulti(): error: other uid %ld not found!\n",other_uid);
		return -1;
	}
	
	if( ! ( (places_map_1i && places_map_1p && places_map_2i && places_map_2p) || only_temporal) ) {
		fprintf(stderr,"crecords::compare_users_onematch_tmulti(): places_maps not set up properly (use the set_compare_params() function before calling this)!\n");
		return -1;
	}
	
	unsigned int mp = 0; // found possible matches
	unsigned int mt = 0; // found temporal matches
	
	uint64_t s2 = it2->second.first;
	uint64_t e2 = it2->second.first + it2->second.second;
	uint64_t s1 = it1->second.first;
	uint64_t e1 = it1->second.first + it1->second.second;
	uint32_t tmaxlast = tmax[ntmax-1];
	
	if(tmin > 0 || tmaxlast > 0) {
		if(tmin > 0) s2 = other_records.tsearch(tmin,true,s2,e2);
		if(s2 < e2 && tmaxlast > 0) e2 = other_records.tsearch(tmaxlast,false,s2,e2) + 1; // + 1, since tsearch returns the result inclusive, but expects the end of the range exclusive
		// adjustment: ts[e-1] == tmaxlast is possible, this needs to be excluded
		while(e2 > s2 && other_records.ts[e2-1] == tmaxlast) e2--; // note: this should not happen if tmaxlast == 0
		
		if(tmin > 0) s1 = tsearch(tmin,true,s1,e1);
		if(s1 < e1 && tmaxlast > 0) e1 = tsearch(tmaxlast,false,s1,e1) + 1;
		// adjustment: ts[e-1] == tmax is possible, this needs to be excluded
		while(e1 > s1 && ts[e1-1] == tmaxlast) e1--;
		
		if(s2 >= e2 || s1 >= e1) {
			for(int i=0;i<ntmax;i++) {
				matches_temporal[i] = 0;
				if(matches_possible) matches_possible[i] = 0;
			}
			return ntmax;
		}
	}
	
	unsigned int dt = dt1;
	if(dt2 > dt) dt = dt2;
	int r = ntmax; // return value (all possible by default, adjusted if impossible matches are found)
	int i = 0; // index in the tmax, matches_temporal and matches_possible arrays
	while(s1 < e1 && s2 < e2) {
		//current records: s1 in this and s2 in other_records
		
		// check that the current records are still in the current time window
		while( ts[s1] >= tmax[i] || other_records.ts[s2] >= tmax[i] ) {
			// save the number of matches found so far
			matches_temporal[i] = mt;
			if(matches_possible) matches_possible[i] = mp;
			
			i++;
			if(i == ntmax) break;
			if(tmax[i] < tmax[i-1] + dt && (tmax[i] > 0 || i < ntmax-1) ) {
				fprintf(stderr,"crecords::compare_users_onematch_tmulti(): tmax sequence is not sorted properly or values too close to each other!\n");
				return -1;
			}
			if(i == r) { // previously we found an impossible match when considering records between [tmax[i-1],tmax[i]) with a record in [tmax[i],tmax[i+1])
				mp = 0;
				only_temporal = true; // no need to check spatial matches anymore
			}
		}
		if(i == ntmax) break;
		
		bool found_impossible = false;
		bool found_possible = false;
		int r1;
		uint32_t impossible_ts = 0; // timestamp of first impossible match found -- it can be < tmax[i] or between tmax[i] and tmax[i+1]
		
		//do the processing for the first one (with respect to temporal order)
		if(ts[s1] < other_records.ts[s2]) {
			if(ts[s1] + dt < other_records.ts[s2]) {
				s1++;
				continue;
			}
			
			r1 = compare_records(s1,other_records,s2,ignore_missing,only_temporal);
			if(r1 == -2) return -1; // error
			if(only_temporal) {
				if(r1) { mt++; s2++; }
				s1++;
				continue;
			}
			else {
				// we have to compare s1 to each record after s2 as well to check for impossible matches
				if(r1 == -1) found_impossible = true;
				else for(uint64_t s22 = s2+1; s22 < e2 && other_records.ts[s22] <= ts[s1] + dt; s22++) {
					int r2 = compare_records(s1,other_records,s22,ignore_missing,false);
					if(r2 == -2) return -1; // error
					if(r2 == -1) {
						found_impossible = true;
						impossible_ts = other_records.ts[s22];
						break;
					}
				}
				
				if(r1 == 1 && (found_impossible == false || (impossible_ts >= tmax[i] && tmax[i] > 0) ) ) {
					found_possible = true;
					for(uint64_t s12 = s1+1; s12 < e1 && ts[s12] <= other_records.ts[s2] + dt; s12++) {
						int r2 = compare_records(s12,other_records,s2,ignore_missing,false);
						if(r2 == -2) return -1; // error
						if(r2 == -1) {
							found_impossible = true;
							if(ts[s12] < impossible_ts || impossible_ts == 0) impossible_ts = ts[s12];
							break;
						}
					}
				}
				if(r1 == -1 || r1 == 1) { s2++; mt++; } // if s1 -- s2 is a temporal match, increase s2 and the temporal match count as well
				s1++;
			}
		}
		else {
			if(ts[s1] > other_records.ts[s2] + dt) {
				s2++;
				continue;
			}
			
			r1 = compare_records(s1,other_records,s2,ignore_missing,only_temporal);
			if(r1 == -2) return -1; // error
			if(only_temporal) {
				if(r1) { mt++; s1++; }
				s2++;
				continue;
			}
			else {
				if(r1 == -1) found_impossible = true;
				else for(uint64_t s12 = s1+1; s12 < e1 && ts[s12] <= other_records.ts[s2] + dt; s12++) {
					int r2 = compare_records(s12,other_records,s2,ignore_missing,false);
					if(r2 == -2) return -1; // error
					if(r2 == -1) {
						found_impossible = true;
						impossible_ts = ts[s12];
						break;
					}
				}
				
				if(r1 == 1 && (found_impossible == false || (impossible_ts >= tmax[i] && tmax[i] > 0) ) ) {
					found_possible = true;
					for(uint64_t s22 = s2+1; s22 < e2 && other_records.ts[s22] <= ts[s1] + dt; s22++) {
						int r2 = compare_records(s1,other_records,s22,ignore_missing,false);
						if(r2 == -2) return -1; // error
						if(r2 == -1) {
							found_impossible = true;
							if(other_records.ts[s22] < impossible_ts || impossible_ts == 0) impossible_ts = other_records.ts[s22];
							break;
						}
					}
				}
				if(r1 == -1 || r1 == 1) { s1++; mt++; } // if s1 -- s2 is a temporal match, increase s1 and the temporal matche count as well
				s2++;
			}
		}
		
		// at this point, mt has been already increased (if needed)
		// what needs to be determined if this is a possible match or impossible
		
		if(found_possible == true && (found_impossible == false || (impossible_ts >= tmax[i] && tmax[i] > 0) ) )
			// this is a possible match and there is no other impossible match (or the impossible match is with a record in the next time window)
			mp++;
		
		if(found_impossible) {
			// check if the impossible match should be counted in this time window or in the next
			if(impossible_ts < tmax[i] || tmax[i] == 0) {
				r = i; // this time window
				mp = 0; // zero out possible matches as wee
				only_temporal = true; // we do not need to consider spatial consistency anymore
			}
			else r = i+1; // next time window (can be the last if there are records beyond tmaxlast
		}
		
	}
	for(;i<ntmax;i++) {
		if(i == r) mp = 0; // if we previously determined that there is an impossible match in this time window, clear possible matches
		matches_temporal[i] = mt;
		if(matches_possible) matches_possible[i] = mp;
	}
	
	return r; // r is the index of the first time window with an impossible match
}






// cluster the records of a user (gived by uid) based on their coordinates
// parameters:
//	tf: optional temporal filter (only include records in the intervals specified there)
//	res: vector to fill with the results; the contents are the central pair and the number of points in that cluster
//	cmin: minimum number of points in a cluster to store in the result
//	rmax: maximum number of clusters to store in the result (zero for these means no limit)
// return value: 0: OK
//	> 0: error (wrong sort order, user was not found)
//	note: if zero points match, or zero cluster are produced that is not reported as an error but results in adding no elements to the result vector
int crecords::cluster_user_coords(int64_t uid, const tfilter* tf, std::vector<std::pair<int,unsigned int> >& res, unsigned int cmin, unsigned int rmax) const {
	if( ! (sorted == uid_cellid || sorted == uid || sorted == uid_time) ) {
		fprintf(stderr,"crecords::cluster_user_coords(): error: records not sorted by user id!\n\t(call sort_by_uid() first)\n");
		return 1;
	}
	
	if(uids_idx_cnts) {
		fprintf(stderr,"crecords::cluster_user_coords(): error: user id index is not set up!\n\t(call sort_by_uid() first)\n");
		return 1;
	}
	
	if(places_map_default == 0) {
		fprintf(stderr,"crecords::cluster_user_coords(): error: no cell neighbor mapping!\n\t(read this from a file first with read_places_map())\n");
		return 1;
	}
	const hashmap<int,hashset<int> >& cell_neighb = *places_map_default;
	
	hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it = uids_idx.find(uid);
	if(it == uids_idx.end()) {
		return 1;
	}
	
	if(it->second.second == 0) {
		return 1;
	}
	
	hashmap<int,unsigned int> cellshist; //create a histogram of the users locations
	
	// fill up the previous hashmap
	for(uint64_t i = it->second.first; i < it->second.first + it->second.second; i++) {
		if(tf) {
			if(tf->filter(ts[i])) cellshist[cells[i]]++;
		}
		else cellshist[cells[i]]++;
	}
	
	unsigned int r = 0;
	
	while(cellshist.size() > 0 && (r < rmax || rmax == 0) ) {
		
		unsigned int cntmax = 0;
		int cellmax = -1;
		for(hashmap<int,unsigned int>::const_iterator it2 = cellshist.begin(); it2 != cellshist.end(); ++it2) {
			int cid = it2->first;
			unsigned int cnt1 = it2->second;
			hashmap<int,hashset<int> >::const_iterator it3 = cell_neighb.find(cid);
			if(it3 != cell_neighb.end()) for(hashset<int>::const_iterator it4 = it3->second.begin(); it4 != it3->second.end(); ++it4) {
				hashmap<int,unsigned int>::const_iterator it5 = cellshist.find(*it4);
				if(it5 != cellshist.end()) cnt1 += it5->second;
			}
			
			if(cnt1 > cntmax) {
				cntmax = cnt1;
				cellmax = cid;
			}
		}
		
		if(cmin) if(cntmax < cmin) break;
		if(cntmax == 0) break;
		
		// add the newly found cluster, remove the cells from the histogram
		res.push_back(std::make_pair(cellmax,cntmax));
		cellshist.erase(cellmax);
		hashmap<int,hashset<int> >::const_iterator it3 = cell_neighb.find(cellmax);
		if(it3 != cell_neighb.end()) for(hashset<int>::const_iterator it4 = it3->second.begin(); it4 != it3->second.end(); ++it4) {
			cellshist.erase(*it4);
		}
	}
	return 0;
}


/*****************************************
 * crecords class, auxilliary functions
 *****************************************/

//functions to ensure that data is properly sorted
		void crecords::sort_by_uid() {
			if( !( sorted == uid || sorted == uid_time || sorted == uid_cellid ) ) {
				if(use_combined_sort) qst3c<int64_t>(&uids,0,(uint64_t)uids.size());
				else {
					if(use_heapsort) heapSort3<int64_t>(&uids,0,(uint64_t)uids.size());
					else qst3<int64_t>(&uids,0,(uint64_t)uids.size());
				}
				sorted = uid;
			}
			cells_idx.clear(); //note: it gets invalidated, make sure we do not have invalid data
			create_uid_idx();
			uids_idx_cnts = false;
		}
		void crecords::sort_by_time() {
			if(sorted != time) {
				if(use_combined_sort) qst3c<uint32_t>(&ts,0,(uint64_t)ts.size());
				else {
					if(use_heapsort) heapSort3<uint32_t>(&ts,0,(uint64_t)uids.size());
					else qst3<uint32_t>(&ts,0,(uint64_t)ts.size());
				}
				sorted = time;
			}
			cells_idx.clear(); //note: it gets invalidated, make sure we do not have invalid data
			uids_idx.clear();
			uids_idx_cnts = false;
		}
		
		void crecords::sort_by_uid_time() {
			if( ! (sorted == uid_time || sorted == uid_time_cellid) ) {
				if(use_combined_sort) qst3c<int64_t,uint32_t>(&uids,0,(uint64_t)ts.size(),&ts);
				else {
					if(use_heapsort) heapSort3<int64_t,uint32_t>(&uids,0,(uint64_t)uids.size(),&ts);
					else qst3<int64_t,uint32_t>(&uids,0,(uint64_t)ts.size(),&ts);
				}
				sorted = uid_time;
			}
			cells_idx.clear(); //note: it gets invalidated, make sure we do not have invalid data
			create_uid_idx();
			uids_idx_cnts = false;
		}
		
		void crecords::sort_by_cellid_time() {
			if( ! (sorted == cellid_time || sorted == cellid_time_uid) ) {
				if(use_combined_sort) qst3c<int,uint32_t>(&cells,0,(uint64_t)cells.size(),&ts);
				else {
					if(use_heapsort) heapSort3<int,uint32_t>(&cells,0,(uint64_t)uids.size(),&ts);
					else qst3<int,uint32_t>(&cells,0,(uint64_t)cells.size(),&ts);
				}
				sorted = cellid_time;
			}
			uids_idx.clear();
			create_cell_idx();
			uids_idx_cnts = false;
		}
		
		void crecords::sort_by_uid_cellid() {
			if(sorted != uid_cellid) {
				if(use_combined_sort) qst3c<int64_t,int>(&uids,0,uids.size(),&cells);
				else {
					if(use_heapsort) heapSort3<int64_t,int>(&uids,0,(uint64_t)uids.size(),&cells);
					else qst3<int64_t,int>(&uids,0,uids.size(),&cells);
				}
				sorted = uid_cellid;
			}
			cells_idx.clear();
			create_uid_idx();
			uids_idx_cnts = false;
		}
		
		void crecords::sort_by_cellid_time_uid() {
			if(sorted != cellid_time_uid) {
				if(use_combined_sort) qst3c<int,uint32_t,int64_t>(&cells,0,nrecords,&ts,&uids);
				else {
					if(use_heapsort) heapSort3<int,uint32_t,int64_t>(&cells,0,nrecords,&ts,&uids);
					else qst3<int,uint32_t,int64_t>(&cells,0,nrecords,&ts,&uids);
				}
				sorted = cellid_time_uid;
			}
			uids_idx.clear();
			uids_idx_cnts = false;
			create_cell_idx();
		}
		
		void crecords::sort_by_uid_time_cellid() {
			if(sorted != uid_time_cellid) {
				if(use_combined_sort) qst3c<int64_t,uint32_t,int>(&uids,0,(uint64_t)ts.size(),&ts,&cells);
				else {
					if(use_heapsort) heapSort3<int64_t,uint32_t,int>(&uids,0,(uint64_t)uids.size(),&ts,&cells);
					else qst3<int64_t,uint32_t,int>(&uids,0,(uint64_t)ts.size(),&ts,&cells);
				}
				sorted = uid_time_cellid;
			}
			cells_idx.clear(); //note: it gets invalidated, make sure we do not have invalid data
			create_uid_idx();
			uids_idx_cnts = false;
		}


//replace cell_id vector according to the given dictionary
		//return: 0: OK, 1: at least one ID was not found
		int crecords::replace_cellids(const hashmap<int,int>& dict) {
			for(size_t i=0;i<cells.size();i++) {
				hashmap<int,int>::const_iterator it = dict.find(cells[i]);
				if(it == dict.end()) {
					fprintf(stderr,"crecords::replace_cellids(): cell id %d was not found in the dictionary supplied!\n",cells[i]);
					return 1;
				}
				cells[i] = it->second;
			}
			if(sorted == cellid_time || sorted == uid_cellid) {
				if(sorted == cellid_time) cells_idx.clear();
				sorted = none;
			}
			return 0;
		}

		
		// return the number of distinct uids
		// note: this recreates the uids_idx hashmap if necessary
		uint64_t crecords::get_uids_count() {
			if(uids.size() == 0) return 0; //no records stored
			if(sorted == uid || sorted == uid_time || sorted == uid_cellid || sorted == uid_time_cellid) {
				if(uids_idx.size() == 0) create_uid_idx();
				return uids_idx.size();
			}
			if(!uids_idx_cnts) uids_count();
			return uids_idx.size();
		}
		
		
		// get a vector with all the user ids among the records
		// the result must be freed with delete by the caller later
		std::vector<int64_t>* crecords::get_uids_p(bool sort) {
			std::vector<int64_t>* r = new std::vector<int64_t>(); //allocate vector structure and initial capacity
			get_uids(r,sort);
			return r;
		}
		// same with returning a (local) vector object (using a move constructor in C++11)
		std::vector<int64_t> crecords::get_uids(bool sort) {
			std::vector<int64_t> r; //allocate vector structure and initial capacity
			get_uids(&r,sort);
			return r;
		}
		// function doing the actual work for the previous two to fill up a user-supplied vector
		void crecords::get_uids(std::vector<int64_t>* r, bool sort) {
			uint64_t cnt = get_uids_count(); //this recreates the hashmap if necessary
			r->clear();
			r->reserve(cnt);
			for(hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it = uids_idx.begin(); it != uids_idx.end(); ++it) {
				r->push_back(it->first);
			}
			if(sort) std::sort(r->begin(), r->end());
		}
		
		// get a vector with the counts of records for each user
		std::vector<std::pair<int64_t,uint32_t> >* crecords::get_uids_cnts_p(bool sort) {
			std::vector<std::pair<int64_t,uint32_t> >* r = new std::vector<std::pair<int64_t,uint32_t> >(); //allocate vector structure and initial capacity
			get_uids_cnts(r,sort);
			return r;
		}
		// same with returning a (local) vector object (using a move constructor in C++11)
		std::vector<std::pair<int64_t,uint32_t> > crecords::get_uids_cnts(bool sort) {
			std::vector<std::pair<int64_t,uint32_t> > r; //allocate vector structure and initial capacity
			get_uids_cnts(&r,sort);
			return r;
		}
		// function doing the actual work for the previous two to fill up a user-supplied vector
		void crecords::get_uids_cnts(std::vector<std::pair<int64_t,uint32_t> >* r, bool sort) {
			uint64_t cnt = get_uids_count(); //this recreates the hashmap if necessary
			r->clear();
			r->reserve(cnt);
			for(hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it = uids_idx.begin(); it != uids_idx.end(); ++it) {
				if(it->second.second > 4294967295UL) // overflow would occur -- note that this would not be undefined behavior
					throw std::runtime_error("crecords::get_uids_cnts(): overflow (at least one user has 2^32 or more records)!\n");
				r->push_back(std::pair<int64_t,uint32_t>(it->first,(uint32_t)(it->second.second)));
					// counts are always stored in the second element
					// use an explicit cast to indicate it should be uint32_t now and also specify the template parameters to make the code more clear
			}
			if(sort) std::sort(r->begin(), r->end());
		}
		
		// get the distribution of number of records per user (i.e. there are r[i].second number of users with r[i].first number of matches)
		std::vector<std::pair<uint32_t,uint32_t> >* crecords::get_records_dist_p() {
			std::vector<std::pair<uint32_t,uint32_t> >* r = new std::vector<std::pair<uint32_t,uint32_t> >();
			get_records_dist(r);
			return r;
		}
		std::vector<std::pair<uint32_t,uint32_t> > crecords::get_records_dist() {
			std::vector<std::pair<uint32_t,uint32_t> > r;
			get_records_dist(&r);
			return r;
		}
		void crecords::get_records_dist(std::vector<std::pair<uint32_t,uint32_t> >* r) {
			uint64_t cnt = get_uids_count(); //this recreates the hashmap if necessary
			hashmap<uint32_t,uint32_t> cnts;
			for(hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it = uids_idx.begin(); it != uids_idx.end(); ++it) {
				if(it->second.second > 4294967295UL) // overflow would occur -- note that this would not be undefined behavior
					throw std::runtime_error("crecords::get_records_dist(): overflow (at least one user has 2^32 or more records)!\n");
				cnts[it->second.second]++; //counts are always stored in the second element, which is implicitly truncated to 32 bits here
			}
			//clear and resize the supplied vector
			r->clear();
			r->reserve(cnts.size());
			//copy results into the supplied vector
			r->insert(r->end(), cnts.cbegin(), cnts.cend());
		}

		// fill the uids_idx hashmap with the number occurrences of each user
		void crecords::uids_count() {
			uids_idx.clear();
			for(uint64_t i=0;i<uids.size();i++) uids_idx[uids[i]].second++; //note: this inserts new elements with the default value (i.e. 0)
			// we store the count in uids_idx.second to be compatible with the version where the uids are actually sorted
			uids_idx_cnts = true;
		}

		// find the position and number of records of one user
		std::pair<uint64_t,uint64_t> crecords::find_user(int64_t uid) const {
			if( ! ( sorted == uid_time || sorted == uid_cellid || sorted == uid_time_cellid || sorted == crecords::uid ) ) {
				fprintf(stderr,"crecords::find_user(): data not sorted by user ID!\n");
				return std::make_pair(nrecords,0);
			}
			
			hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it1 = uids_idx.find(uid);
			if(it1 == uids_idx.end()) {
				fprintf(stderr,"crecords::find_user(): error: uid %ld not found!\n",uid);
				return std::make_pair(nrecords,0);
			}
			
			return it1->second;
		}
		
		
		// get position or iterator for a user (similar to the previous one, but slightly different interface)
		// get the first record with the given uid (or nrecords if not found)
		uint64_t crecords::get_uid_first(int64_t uid) const {
			if( ! (sorted == uid || sorted == uid_cellid || sorted == uid_time || sorted == uid_time_cellid ) ) {
				fprintf(stderr,"crecords::get_uid_first(): records not sorted by uid!\n");
				return nrecords;
			}
			if(uids_idx_cnts || (uids_idx.size() == 0 && uids.size() > 0) ) {
				fprintf(stderr,"crecords::get_uid_first(): no index for user IDs!\n");
				return nrecords;
			}
			hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it = uids_idx.find(uid);
			if(it == uids_idx.end()) return nrecords;
			else return it->second.first;
		}
		
		// get the last+1 record with the given uid (or nrecords if not found -- i.e. that is not necessarily an error, but the return value
		//		of get_uids_first() should be checked for that first)
		uint64_t crecords::get_uid_last(int64_t uid) const {
			if( ! (sorted == uid || sorted == uid_cellid || sorted == uid_time || sorted == uid_time_cellid ) ) {
				fprintf(stderr,"crecords::get_uid_last(): records not sorted by uid!\n");
				return nrecords;
			}
			if(uids_idx_cnts || (uids_idx.size() == 0 && uids.size() > 0) ) {
				fprintf(stderr,"crecords::get_uid_last(): no index for user IDs!\n");
				return nrecords;
			}
			hashmap<int64_t,std::pair<uint64_t,uint64_t> >::const_iterator it = uids_idx.find(uid);
			if(it == uids_idx.end()) return nrecords;
			else return it->second.first + it->second.second;
		}
		
		// get an iterator starting at the records of the given user ID
		crecords::const_iterator crecords::get_uid_iterator_begin(int64_t uid) const {
			uint64_t i = get_uid_first(uid);
			if(i >= nrecords) {
				fprintf(stderr,"crecords::get_uid_iterator_begin(): user id %ld not found or error searching!\n",uid);
				throw new std::runtime_error("crecords::get_uid_iterator_begin(): error searching!\n");
			}
			return crecords::const_iterator(this,i);
		}
		
		// get a sentinel for the last record of a given user ID
		crecords::sentinel crecords::get_uid_iterator_end(int64_t uid) const {
			uint64_t i = get_uid_last(uid);
			return crecords::sentinel(i); // no need to check for error, get_uid_iterator_begin already throws an exception for that
		}

// "implementation" of getting iterators out of a crecords class
crecords::const_iterator crecords::begin() const {
	return crecords::const_iterator(this);
}
crecords::const_iterator crecords::cbegin() const {
	return crecords::const_iterator(this); // all iterators are constant iterators, the records can only be changed with the add/load functions
}

crecords::sentinel crecords::end() const { // this should be used instead of an end iterator; it will not work with C++ std algorithms but will work in a for loop
	return crecords::sentinel(nrecords);
}
crecords::sentinel crecords::cend() const {
	return crecords::sentinel(nrecords);
}	


#ifdef USECEREAL
		//read data from a specified file
		int crecords::read_crecords_serialized(crecords& cr, const char* fn, bool zip) {
			if(zip) {
				struct popen_noshell_pass_to_pclose spc;
				const char* const popen_args[] = { "/bin/gzip", "-cd", fn, 0 };
				FILE* f = popen_noshell(popen_args[0],popen_args,"r",&spc,0);
				if(f == 0) return 1;
				stdiostream sifs(f);
				cereal::PortableBinaryInputArchive ar(sifs.stream());
				ar(cr);
				pclose_noshell(&spc);
			}
			else {
				std::ifstream ifs(fn,std::ifstream::binary);
				cereal::PortableBinaryInputArchive ar(ifs);
				ar(cr);
			}
			return 0;
		}
		
		//write data to a specified file
		int crecords::write_crecords_serialized(const crecords& cr, const char* fn) {
			//save cr using the serialization
			std::ofstream of(fn,std::ofstream::trunc | std::ofstream::binary);
			cereal::PortableBinaryOutputArchive ar(of);
			ar(cr);
			return 0;
		}
#endif


// for debugging
int crecords::places_map_count(const hashmap<int,hashset<int> >& places_map, int place_id) {
	return places_map.count(place_id);
}
int crecords::places_map_count(const hashmap<int,hashset<int> >& places_map, int place_id, int cell_id) {
	hashmap<int,hashset<int> >::const_iterator it = places_map.find(place_id);
	if(it == places_map.end()) return 0;
	else return it->second.count(cell_id);
}

		

/**********************************************************
 * separate class for storing the temporal distribution
 * of points from one of the datasets
 * can be generated from the crecords class
 **********************************************************/

// generate random timestamps based on the distributions
// the given seed is used to initialize the random number generator if useseed == true
//	(otherwise, the default initialization in the implementation is used)
// return 0 on success, 1 on error (memory allocation failure or if the PDFs were not initialized previously)
int tdist::rg_prepare(bool useseed, uint64_t seed) {
	if( tmax < tmin ) return 1;
	scdf.clear();
	
	scdf.resize(tmax-tmin+1);
	
	uint64_t ssum = 0;
	for(uint32_t i=0; i <= (tmax-tmin); i++) {
		ssum += sdist[i];
		scdf[i] = ssum;
	}
	
	if(useseed) rg.seed(seed);
	
	return 0;
}

// get a random timestamp according to the distribution from the given range (max is inclusive)
// if min == max == 0, the whole time range in this instance is used
// if the range is invalid or the cdfs are uninitialized and exception is thrown
// if no timestamps are in the range, a timestamp out of the range is returned (0, if it is not in the range)
uint32_t tdist::get_random_ts_r(std::mt19937_64& rg1, uint32_t min, uint32_t max) const {
	if( scdf.size() < (tmax-tmin+1) ) throw new std::runtime_error("tdist::get_random_ts(): CDF arrays not prepared (use tdist::rg_prepare() first)!\n");
	if(max < min) throw new std::runtime_error("tdist::get_random_ts(): invalid parameters!\n");
	if(min == 0 && max == 0) {
		min = tmin;
		max = tmax;
	}
	if(min == max) return min;
	
	// generate a random number between scdf[min-1] and scdf[max], the result will be the smallest t such that scdf[t] >= r
	uint64_t rmin = 0;
	if(min > tmin) rmin = scdf[min-tmin-1]; // else 0
	uint64_t rmax = scdf[max-tmin-1]; //note: max >= tmin+1 here, otherwise max == min, which was taken care of previously
	
	uint64_t r1 = rmin + get_random_range(rg1,rmax-rmin);
	uint32_t r2 = std::lower_bound(scdf.begin(), scdf.end(), r1) - scdf.begin();
	
	return r2+tmin;
}


/********************************************************
 * separate class for storing user record distributions	
 * and generating user IDs randomly who have the given
 * number of records
 ********************************************************/	
void user_rand::sort_uids() {
	std::sort(uids.begin(), uids.end(), [](const std::pair<int64_t,uint64_t> &left, const std::pair<int64_t,uint64_t> &right) {
	    if(left.second < right.second) return true;
	    if(left.second > right.second) return false;
	    return left.first < right.first;
	});
	nusers = uids.size();
	if(nusers > 0) { cntmax = uids[nusers-1].second; cntmin = uids[0].second; }
	else { cntmax = 0; cntmin = 0; }
}
		

int64_t user_rand::get_random_uid(uint32_t min, uint32_t max, uint32_t* nrecords) {
	if(uids.size() == 0) return -1;
	if(min == 0 && max == 0) { min = cntmin; max = cntmax; }
	if(min > cntmax || max < cntmin) return -1;
	// find min and max using binary search in C++ STL
	uint32_t imin2 = 0;
	uint32_t imax2 = uids.size();
	
	std::vector<std::pair<int64_t, uint32_t> >::const_iterator imin,imax;
	if(min > cntmin) { // minor optimization: avoid search if we know that the first or last element is needed
		imin = std::lower_bound(uids.begin(),uids.end(),std::pair<int64_t,uint32_t>(-1,min),
			[](const std::pair<int64_t,uint32_t> &left, const std::pair<int64_t,uint32_t> &right) { return left.second < right.second; });
		imin2 = imin - uids.cbegin(); // easier to deal with numbers than iterators
	}
	if(max < cntmax) {
		imax = std::upper_bound(uids.begin(),uids.end(),std::pair<int64_t,uint32_t>(-1,max),
			[](const std::pair<int64_t,uint32_t> &left, const std::pair<int64_t,uint32_t> &right) { return left.second < right.second; }); 
		imax2 = imax - uids.cbegin(); // easier to deal with numbers than iterators
	}
	if(imin2 >= imax2 || imin2 >= uids.size()) return -1; //no user in the given range
	
	// choose a user randomly from imax2 - imin2 users
	uint32_t i = imin2 + get_random_range(r,imax2 - imin2);
	if(nrecords) *nrecords = uids[i].second;
	return uids[i].first;
}


// initialize the length distributions (read from the given file)
int tlendistr::init_length_distr(FILE* f) {
	// we expect a length and the corresponding frequency on each line
	// the lengths should be sorted
	lengths.clear();
	lcdf.clear();
	
	unsigned int line = 0;
	uint32_t len0 = 0;
	uint64_t sum = 0;
	while(1) {
		line++;
		uint32_t len;
		uint64_t freq;
		int a = fscanf(f,"%u %lu",&len,&freq); // note: newlines are skipped, so this might accept bad format (one number per line) as well
		if(a == EOF) break;
		if(a != 2) {
			fprintf(stderr,"tdistrpairiterator::init_length_distr(): invalid data on line %u of input!\n",line);
			return 1;
		}
		if(len <= len0 && !(len == 0 && len0 == 0 && lengths.size() == 0)) {
			fprintf(stderr,"tdistrpairiterator::init_length_distr(): input not sorted on line %u!\n",line);
			return 1;
		}
		sum += freq;
		lengths.push_back(len);
		lcdf.push_back(sum);
		do {
			a = getc(f);
		} while( ! ( a == '\n' || a == EOF ) );
		if(a == EOF) break;
	}
	return 0;
}




