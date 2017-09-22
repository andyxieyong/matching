/*
 * cdrrecords.h
 * 
 * read CDR or similar checkin records, perform spatiotemporal search, various helpers
 * 
 * note: this requires a C++11 compiler (tested with g++ 4.8, 5.3, clang 3.4)
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


#ifndef CDRRECORDS_H
#define CDRRECORDS_H

#define __STDC_LIMIT_MACROS

#include <stdio.h> //std file io
#include <stdint.h> //for uint64_t definition
#include <vector> //storing dynamic data
#include <deque> //used by the temporal distribution class
#include <utility> //std::pair
#include <algorithm> //std::sort, std::lower_bound, std::upper_bound


#ifdef USECEREAL //include support for the cereal c++ serialization library -- this implies using std::unordered_map containers
#include <fstream>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/unordered_set.hpp>
#include <cereal/types/utility.hpp>
// for reading gzipped binary input (by creating a gzip process with popen_noshell())
#include "popen_noshell.h"
#endif

// use hashmaps and hashsets provided by the C++11 STL
// define a shorthand type name for these
// (this allows for shorter code and switching to an alternative compatible implementation)
#include <unordered_map>
#include <unordered_set>
template <class T, class U> using hashmap = std::unordered_map<T,U>;
template <class T> using hashset = std::unordered_set<T>;

// random generator
// #include "mt19937.h"
#include <random>


/**************************
 * auxilliary data types
 **************************/


// store one record in one struct -- note: this is not used internally, only for interfacing with callers
struct record {
	int64_t uid;
	uint32_t ts;
	int cell_id;
	int8_t startstop; // transit start or stop: 0: start, 1: end, -1: not applicable
};


// temporal filter definition
struct tfilter {
	unsigned int start; //filter start (seconds in the day)
	unsigned int end; //filter end (i.e. times between start and end are considered each day)
	bool weekend; // if false, all times from weekend days are omitted
	bool weekendall; // if true (and weekend is true also), all times from weekend days are included
	bool invert; // if true, filter interval is inverted (but not the weekend switch, i.e. to include weekday and weekend nights, use weekend == true
		// and invert == true)
	
	//determine if a timestamp falls in the times as specified by this filter
	inline bool filter(time_t t) const {
		struct tm t0;
		struct tm* t1 = gmtime_r(&t,&t0);
		if(!t1) return false;
		
		if(!weekend) if(t1->tm_wday == 0 || t1->tm_wday == 6) return false;
		if(weekendall) if(t1->tm_wday == 0 || t1->tm_wday == 6) return true;
		
		unsigned int seconds = 3600 * t1->tm_hour + 60 * t1->tm_min + t1->tm_sec;
		
		if(invert) {
			if(seconds < start || seconds >= end) return true;
			else return false;
		}
		else {
			if(seconds < start || seconds >= end) return false;
			else return true;
		}
	}
};


// previous functions for temporal filtering for night / day; times are considered hours of day here

//determine if a timestamp is night (1) or day (2) or neither (0)
inline int dnt(time_t t, unsigned int daystart, unsigned int dayend, unsigned int nightstart, unsigned int nightend) {
	struct tm* t1 = gmtime(&t);
	if(!t1) return 0;
	if(t1->tm_hour < dayend && t1->tm_hour >= daystart) return 2;
	if(t1->tm_hour < nightend || t1->tm_hour >= nightstart) return 1;
	return 0;
}

//determine if a timestamp is night (1) or day (2), weekend night (3) or day (4) or neither (0) or weekend in general (5)
inline int dnt2(time_t t, unsigned int daystart, unsigned int dayend, unsigned int nightstart, unsigned int nightend) {
	struct tm* t1 = gmtime(&t);
	if(!t1) return 0;
	if(t1->tm_hour < dayend && t1->tm_hour >= daystart) {
		if(t1->tm_wday == 0 || t1->tm_wday == 6) return 4;
		else return 2;
	}
	if(t1->tm_hour < nightend || t1->tm_hour >= nightstart) {
		if(t1->tm_wday == 0 || t1->tm_wday == 6) return 3;
		else return 1;
	}
	if(t1->tm_wday == 0 || t1->tm_wday == 6) return 5;
	return 0;
}



/*******************************************************************************
 * main class for storing "records" of users (either CDR or transportation
 * 	or possibly some other dataset with a finite number of placed of interest)
 *******************************************************************************/

class crecords_iterator; // forward declaration of iterator and sentinel to be used by the class
class crecords_iterator_sentinel;
class crecords {
	public:
		enum sortorder { none, uid, time, uid_time, cellid_time, uid_cellid, cellid_time_uid, uid_time_cellid };
		enum places_map_type { place_cell, place_place, cell_cell, cell_cell_centers, cell_place };
		
	protected:
	//main data: user_id (/number), timestamp, cell id
		std::vector<int64_t> uids;
		std::vector<uint32_t> ts;
		std::vector<int> cells;
		std::vector<bool> startstop; //startstop[i] == false if the ith record is the start of a trip, true if it is the end (makes most sense if
			// records are ordered by user ID and time; this is not necessary required to be stored / utilized (only for the transportation records)
		uint64_t nrecords; //total number of records stored in the previous arrays (should be equal to uids.size(), ts.size() and cells.size()!)
		
		// store place (stop) <-> cell matchings here -- using several possible radii
		std::vector<hashmap<int,hashset<int> > > places_maps;
		std::vector<double> radii;
		places_map_type pmt; // type of data saved in places_map
		
		//helper to store indices in cells to help efficient search
		hashmap<int,std::pair<uint64_t,uint64_t> > cells_idx;
		//same for user ids (if records are sorted by uid)
		hashmap<int64_t,std::pair<uint64_t,uint64_t> > uids_idx;
		bool uids_idx_cnts; //if true, the uids_idx hashmap is filled up, but contains the number of records per user and not actual indices
			// (this way, the hash map can be used to count the number of users)
		// otherwise the hashmap is either empty or contains indices
		
		
		sortorder sorted;
		
		// double radius; //radius used to create the above mapping -- store this along with the map, so it can be better reused
		
		
		// search parameters
		const hashmap<int,hashset<int> >* places_map_default; // default map to use for searching (for the search_records*() functions)
		
		// compare_users*() functions:
		unsigned int dt1; // time period for "walking": used before a record if startstop == false, used after a record when startstop == true
		unsigned int dt2; // time period for "transit": used before a record if startstop == true, used after a record when startstop == false
		const hashmap<int,hashset<int> >* places_map_1p; // places reachable with "walking": 1p contains places which are considered possible,
		const hashmap<int,hashset<int> >* places_map_1i; // while 1i contains plaves which are impossible
		const hashmap<int,hashset<int> >* places_map_2p; // places reachable with "transit": similarly to the previous
		const hashmap<int,hashset<int> >* places_map_2i; // these are used the same way as dt1 and dt2
		
		
	public:
		crecords() {
			nrecords = 0;
			sorted = none;
			uids_idx_cnts = false;
			use_heapsort = false;
			use_combined_sort = true;
			places_map_default = 0;
			timeoverlap = true;
			
			places_map_1p = 0;
			places_map_1i = 0;
			places_map_2p = 0;
			places_map_2i = 0;
			dt1 = 0;
			dt2 = 0;
		}
		
		
		/****************************************************************************************
		 * 1. functions for reading / loading data
		 ****************************************************************************************/
		// read records from TSV data, format should be uid\tts\tcell_id with possibly more stuff after that
		// return value: number of records read
		// note: it can throw an exception if it runs out of memory (or a SIGSEGV / OOM killer on Linux if memory overcommiting is enabled)
		//~ uint64_t records_read_tsv(FILE* in, unsigned int header_skip = 1, uint64_t lines_max = 0);
		
		// load / add records from a custom data source, supplied by the caller (the load_ version clears all data, the add_ version extends)
		// the iterator supplied as it should be dereferenceable as const record&, and incrementable
		// the sentinel needs to be comparable to the iterator and return true if the end is reached (i.e. it can be the same class as iterator
		//		or can be a "sentinel" class which checks if the iterator is still valid in line with the range concept:
		//		http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4128.html or C# IEnumerables)
		// records are loaded until it == end
		// if nrecords > 0, at maximum this many records are loaded
		// if ignore_end == true, the end iterator is not accessed (no comparisons performed), exactly nrecords records are loaded
		//		(useful for iterators which do not traverse a finite set, e.g. generate random data)
		// iterators which can be useful are provided later, e.g. tsv_iterator, crecords_iterator, filter_iterator, randomiterator and tdistiterator
		//		these can be used for loading data with an optional filtering or random data
		//		see sertest.cpp and random_match.cpp for usage examples
		template<class iterator, class sentinel> uint64_t add_records_custom(iterator& it, const sentinel& end, uint64_t nrecords = 0, bool ignore_end = false) {
			bool sorted_by_time = false;
			int8_t sorted_by_uid = 0; // 3: uid, time and cellid; 2: uid and time; 1: only uid; 0: not sorted
			bool sorted_by_uid_cid = false;
			int8_t sorted_by_cid = 0; // 3: cellid, time and uid; 2: cellid and time; 1: only cellid; 0: not sorted
			if(this->nrecords > 0) switch(sorted) {
					// possible sort orders: { none, uid, time, uid_time, cellid_time, uid_cellid, cellid_time_uid, uid_time_cellid };
				case time:	sorted_by_time = true; break;
				case uid:	sorted_by_uid = 1; break;
				case uid_time:	sorted_by_uid = 2; break;
				case cellid_time: sorted_by_cid = 2; break;
				case uid_cellid: sorted_by_uid_cid = true; break;
				case cellid_time_uid: sorted_by_cid = 3; break;
				case uid_time_cellid: sorted_by_uid = 3; break;
			}
			else {
				sorted_by_time = true;
				sorted_by_uid = 3;
				sorted_by_cid = 3;
				sorted_by_uid_cid = true;
			}
			
			uids_idx.clear();
			uids_idx_cnts = false;
			cells_idx.clear();
			
			if(nrecords) reserve(this->nrecords + nrecords);
			uint64_t i = 0;
			for(;;++it,i++) {
				if(ignore_end) { if(i >= nrecords) break; } // note: these conditions seemed too complex / long to include in the for() statement
				else { if( it == end || (nrecords > 0 && i >= nrecords) ) break; }
				record r = *it;
				
				if(r.startstop == 0 || r.startstop == 1) {
					if(startstop.size() == 0 && uids.size() > 1) {
						fprintf(stderr,"crecords::add_records_custom(): start / stop value supplied at record %lu while not previously!\n",i);
						return i;
					}
					if(r.startstop == 0) startstop.push_back(false);
					else startstop.push_back(true);
				}
				else if(startstop.size() > 0) {
					fprintf(stderr,"crecords::add_records_custom(): start / stop value missing at record %lu!\n",i);
					return i;
				}
				
				// test if data is still sorted
				if(i > 0 || this->nrecords > 0) {
					if(sorted_by_time) if(r.ts < ts.back()) sorted_by_time = false;
					if(sorted_by_uid_cid) if(r.uid < uids.back() || ( r.uid == uids.back() && r.cell_id < cells.back() ) ) sorted_by_uid_cid = false;
					if(sorted_by_uid) {
						if(r.uid < uids.back()) sorted_by_uid = 0;
						else if(sorted_by_uid > 1 && r.uid == uids.back()) {
							if(r.ts < ts.back()) sorted_by_uid = 1;
							else if(sorted_by_uid > 2 && r.ts == ts.back()) if(r.cell_id < cells.back()) sorted_by_uid = 2;
						}
					}
					if(sorted_by_cid) {
						if(r.cell_id < cells.back()) sorted_by_cid = 0;
						else if(sorted_by_cid > 1 && r.cell_id == cells.back()) {
							if(r.ts < ts.back()) sorted_by_cid = 1;
							else if(sorted_by_cid > 2 && r.ts == ts.back()) if(r.uid < uids.back()) sorted_by_cid = 2;
						}
					}
				}
				
				uids.push_back(r.uid);
				ts.push_back(r.ts);
				cells.push_back(r.cell_id);
			}
			this->nrecords = uids.size();
			if(i > 0) {
				// if records were added, check if sort order needs to be adjusted
				sorted = none;
				if(sorted_by_time) sorted = time;
				if(sorted_by_uid_cid) sorted = uid_cellid;
				switch(sorted_by_uid) {
					case 1: sorted = uid; break;
					case 2: sorted = uid_time; break;
					case 3: sorted = uid_time_cellid; break;
				}
				switch(sorted_by_cid) {
					case 1: sorted = none; // not recorded (not useful)
					case 2: sorted = cellid_time;
					case 3: sorted = cellid_time_uid;
				}
			}
			return i;
		}
		template<class iterator, class sentinel> uint64_t load_records_custom(iterator& it, const sentinel& end, uint64_t nrecords = 0, bool ignore_end = false) {
			uids.clear();
			ts.clear();
			cells.clear();
			startstop.clear();
			return add_records_custom(it,end,nrecords,ignore_end);
		}
		
		
		// auxilliary function to reserve storage in advance (useful if autogrowth would result in too much memory allocation)
		void reserve(uint64_t size) { uids.reserve(size); ts.reserve(size); cells.reserve(size); }
		
		// clear all data stored (note: this might not actually free memory due to the way std::vector works; the only sure way to free memeory is
		//		to destruct the class)
		void clear() { uids.clear(); ts.clear(); cells.clear(); uids_idx.clear(); cells_idx.clear(); startstop.clear(); uids_idx_cnts = false; sorted = none; }
		
		// shrink space used by the vectors to store data (note: this might not actually free memory due to the way std::vector works; the only
		//		sure way to free memeory is to destruct the class)
		void shrink() { uids.shrink_to_fit(); ts.shrink_to_fit(); cells.shrink_to_fit(); startstop.shrink_to_fit(); }
		
		
		/******************************************************************
		 * 2. functions for searching among the records or counting the
		 * 	number of matches between a given pair of users
		 ******************************************************************/
		// search for events among the records
		// process all events in the timespan [tmin,tmax] (inclusive), add possible / impossible matches to the containers supplied
		//		(note: the containers are not emptied before, that is the responsibility of the caller)
		//		the containers need to implement a method end() and insert(), to be used as c.insert(c.end(), int64_t id) (vector and unordered_set do this)
		// method: process all records in the timespan, all of them (i.e. the user ID from them) go into either to the possible / impossible
		//		(note: this could later be improved by a "buffer zone")
		// decision is made based on places_map: for each event it is considered "possible" if it is in places_map[place_id] and
		//		impossible otherwise (note: it is an error if place_id is not found!)
		// note: a uid can end up in both sets (if it appears in more than one event in the time span), then the caller can decide
		//		how to further process it
		// if only_temporal == true, spatial consistency is not checked and all matches are added to possible_matches
		// note: the two containers can be the same (if possible / impossible is not important, in this case setting also only_temporal == true will
		//		speed up processing as well)
		// return value: 0 if search was OK, 1 if an error occured (place_id not found, or tmin > tmax) -- note: I'm not familiar with C++ exceptions,
		//		but we could just use them in this case too
		// suggestion: maybe return the number of total records processed? that should be evident from the vectors' size, too
		// note: this expects the records to be in the correct sorted order and return an error if this is not the case 
		
		template<class container>
		int search_records(int place_id, unsigned int tmin, unsigned int tmax, container& possible_matches, container& impossible_matches,
				bool only_temporal = false) const {
			if(uids.size() == 0) return 1;
			if(tmin > tmax) return 1;
			if(sorted != time) {
				fprintf(stderr,"crecords::search_records(): error: data not sorted by time!\n");
				return 1;
			}
			if(places_map_default == 0 && only_temporal == false) {
				fprintf(stderr,"crecords::search_records(): error: places_map not set!\n");
				return 1;
			}
			uint64_t start = tsearch(tmin,true);
			if(start == nrecords) return 0; //no records in the time range
			uint64_t end = tsearch(tmax,false);
			if(end == nrecords) return 0; //no records in the time range (note: this can occur if tmin <= tmax < ts[0], in this case start == 0)
			
			if(only_temporal) for(uint64_t i=start;i<=end;i++)
				possible_matches.insert(possible_matches.end(), uids[i]);
			else {
				const hashmap<int,hashset<int> >& places_map = *places_map_default;
				hashmap<int,hashset<int> >::const_iterator it = places_map.find(place_id);
				if(it == places_map.end()) return 1;
				const hashset<int>& cellmatches = it->second;
				
				for(uint64_t i=start;i<=end;i++) {
					int cell = cells[i];
					if(cellmatches.count(cell) == 0) impossible_matches.insert(impossible_matches.end(), uids[i]);
					else possible_matches.insert(possible_matches.end(), uids[i]);
				}
			}
			return 0;
		}
		
		
		// search records but store only possible matches (similarly to the previous one)
		// use records sorted by cell_id and timestamp
		// note: this expects the records to be in the correct sorted order and return an error if this is not the case 
		template<class container>
		int search_records_possible(int place_id, unsigned int tmin, unsigned int tmax, container& possible_matches) const {
			if(uids.size() == 0) return 1;
			if(tmin > tmax) return 1;
			if( ! (sorted == cellid_time || sorted == cellid_time_uid) ) {
				fprintf(stderr,"crecords::search_records_possible(): error: data not sorted by cell_id,time!\n");
				return 1;
			}
			if(places_map_default == 0) {
				fprintf(stderr,"crecords::search_records(): error: places_map not set!\n");
				return 1;
			}
			const hashmap<int,hashset<int> >& places_map = *places_map_default;
			
			hashmap<int,hashset<int> >::const_iterator it = places_map.find(place_id);
			if(it == places_map.end()) return 1;
			const hashset<int>& cellmatches = it->second;
			for(hashset<int>::const_iterator it2 = cellmatches.begin(); it2 != cellmatches.end(); ++it2) {
				hashmap<int,std::pair<uint64_t,uint64_t> >::const_iterator it3 = cells_idx.find(*it2);
				if(it3 != cells_idx.end()) {
					uint64_t s = it3->second.first;
					uint64_t e = s + it3->second.second;
					
					uint64_t start = tsearch(tmin,true,s,e);
					if(start >= e) continue; //not found (no CDR events in this cell in this time window)
					uint64_t end = tsearch(tmax,false,s,e);
					if(end >= e) continue;
					
					//add all uids between start and end (inclusive) to the set of possible matches
					for(uint64_t i=start;i<=end;i++) possible_matches.insert(possible_matches.end(), uids[i]);
				}
			}
			return 0;
		}
		
		
		// compare uid in this dataset with other_uid in the dataset other_records
		// a possible match is considered if for a record i of uid (i.e. uids[i], cells[i], ts[i]), if other_uid has a 
		//	record j in other_records (i.e. other_records.uid[j], .cells[j], .ts[j]), where ts[i] - dt1 <= other_recods.ts[j] <= ts[i] + dt2
		//	and other_records.places_map[ cells[i] ].count( other_records.cells[j] ) > 0
		// search is limited in the range [tmin,tmax]
		// an impossible match is considered if the previous count() is 0
		// if ignore_missing == true, a missing place ID results in an impossible match, not an error
		// return value: 1: possible (*matches is the number of matches), 0: impossible (at least one inconsistent match), -1: error while searching (place_id not found)
		// note: this expects the records to be in the correct sorted order (in both crecords classes) and returns an error if this is not the case 
		inline int compare_users(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches,
			bool ignore_missing = false, uint32_t tmin = 0, uint32_t tmax = 0) const {
				return compare_users2(uid,other_records,other_uid,matches,0,false,ignore_missing,false,tmin,tmax);
		}
		
		// similar to the previous one, but each record (in each dataset) can be matched only once (at maximum)
		// thus the number of matches is at maximum the smaller of the number of records the two users have
		// this is achieved by processing records of the two users in temporal order and in parallel and when multiple possible records
		//	match, selecting the one with the minimum timestamp as the "real" match
		// search is limited in the range [tmin,tmax]
		// return value: 1: possible (*matches is the number of matches), 0: impossible (at least one inconsistent match), -1: error while searching (place_id not found)
		inline int compare_users_onematch(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches,
			bool ignore_missing = false, uint32_t tmin = 0, uint32_t tmax = 0) const {
				return compare_users2_onematch(uid, other_records, other_uid, matches, 0, false, ignore_missing, false,tmin,tmax);
		}
			
		// compare two users but do not consider spatial consistency, i.e. count only matches in time without considering spatial consistency
		// return value: 1: OK, -1: one of the users was not found
		inline int compare_users_temporal(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches,
			uint32_t tmin = 0, uint32_t tmax = 0) const {
				return compare_users2(uid, other_records, other_uid, matches, 0, true, false, false,tmin,tmax);
		}
		// compare two users but do not consider spatial consistency, i.e. count only matches in time without considering spatial consistency
		// only consider records from each dataset only once
		// return value is the number of matches
		inline int compare_users_temporal_onematch(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches,
			uint32_t tmin = 0, uint32_t tmax = 0) const {
				return compare_users2_onematch(uid, other_records, other_uid, matches, 0, true, false, false,tmin,tmax);
		}
		
		// compare two users and return the number of both possible and impossible matches
		// return value: 1: possible, 0: impossible (at least one inconsistent match), -1: error while searching
		inline int compare_users_all(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches_possible,
			unsigned int* matches_impossible, bool ignore_missing = false, uint32_t tmin = 0, uint32_t tmax = 0) const {
				return compare_users2(uid, other_records, other_uid, matches_possible, matches_impossible, false, ignore_missing, true, tmin, tmax);
		}
		// compare two users and return the number of both possible and impossible matches
		// count each event in both datasets only once -- note: if an event would be counted as possible, but also has a temporal match
		//	which is spatially impossible, this might be problematic !!
		// return value: 1: possible, 0: impossible (at least one inconsistent match), -1: error while searching
		inline int compare_users_all_onematch(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches_possible,
			unsigned int* matches_impossible, bool ignore_missing = false, uint32_t tmin = 0, uint32_t tmax = 0) const {
				return compare_users2_onematch(uid, other_records, other_uid, matches_possible, matches_impossible,
					false, ignore_missing, true,tmin,tmax);
		}
		
		
		// functions doing the actual work for all the above
		// two versions: one for counting all, one for considering each event only once
		
		// note: these functions now only use the places_maps (and the time paraneters) from this class as set up by set_compare_params()
		// if trip start / stop information is present, then these functions should be called from the class containing them
		
		// count all matches for events
		// number of possible matches returned in *matches_possible
		// if only_temporal is true, spatial consistency is not checked and *matches_possible contains all temporal matches
		// if count_impossible == false (default), than one impossible match aborts the compaision and 0 is returned (i.e. matches_possible and
		//		matches_impossible are not updated); if all matches are possible, then *matches_possible contains the number of possible (temporal)
		//		matches, but matches_impossible is not dereferenced (i.e. is count_impossible == false, then matches_impossible can be NULL or invalid)
		// if count_impossible == true, and there is at least one impossible match then *matches_possible contains the number of temporal matches (as if
		//		only_temporal == true was given); if matches_impossible is not NULL, than the number of impossible matches is counted and stored there
		// return value: 1: possible (*matches is the number of matches), 0: impossible (at least one inconsistent match), -1: error while searching (place_id not found)
		int compare_users2(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches_possible,
			unsigned int* matches_impossible, bool only_temporal, bool ignore_missing, bool count_impossible, uint32_t tmin = 0, uint32_t tmax = 0) const;
		// count matches, but count each event from each dataset only once maximum
		// number of possible matches returned in *matches_possible
		// if only_temporal is true, spatial consistency is not checked and *matches_possible contains all temporal matches
		// if count_impossible == false (default), than one impossible match aborts the compaision and 0 is returned (i.e. matches_possible and
		//		matches_impossible are not updated); if all matches are possible, then *matches_possible contains the number of possible (temporal)
		//		matches, but matches_impossible is not dereferenced (i.e. is count_impossible == false, then matches_impossible can be NULL or invalid)
		// if count_impossible == true, and there is at least one impossible match then *matches_possible contains the number of temporal matches (as if
		//		only_temporal == true was given); if matches_impossible is not NULL, than the number of impossible matches is counted and stored there
		// return value: 1: possible (*matches is the number of matches), 0: impossible (at least one inconsistent match), -1: error while searching (place_id not found)
		int compare_users2_onematch(int64_t uid, const crecords& other_records, int64_t other_uid, unsigned int* matches_possible,
			unsigned int* matches_impossible, bool only_temporal, bool ignore_missing, bool count_impossible, uint32_t tmin = 0, uint32_t tmax = 0) const;
		
		
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
		int compare_users_onematch_tmulti(int64_t uid, const crecords& other_records, int64_t other_uid, int ntmax,
			unsigned int* matches_temporal, unsigned int* matches_possible, uint32_t tmin, const uint32_t* tmax,
			bool only_temporal, bool ignore_missing) const;
		
		// set the parameters used for comparisons by the above functions
		// need to be called before calling any of those, but after the coordinates and the Voronoi-tesselation has been loaded
		// not thread-safe (since it modifies the stored params), should be called by the main thread at first
		// generate appropriate places_maps if needed
		// if searchreverse == true, the mappings are created reversed (cell_id -> place_id), this should be used if this class contains transportation records
		int set_compare_params(unsigned int dt1_, unsigned int dt2_, double radius1p, double radius1i, double radius2p, double radius2i,
			bool searchreverse = false);
		
		// further parameter: can timespans for recrods overlap? (true by default)
		bool timeoverlap;
		
	protected:
		// helper function for the above
		int compare_records(uint64_t i, const crecords& other_records, uint64_t other_i, bool ignore_missing, bool only_temporal) const;
	
	public:
		
		/***********************************************************
		 * 3. functions for reading place-cell mappints
		 * (note: generating these has been moved to voronoi_map_create.cpp)
		 ***********************************************************/
		
		
		inline places_map_type get_places_map_type() { return pmt; }
		
		// read a place-cell mapping from the specified file and add it to the list (places_maps)
		// return values: 0: success; >0: error
		// it is considered an error if the type is different than what was previously read or if the radius is exactly the same as previously
		int read_places_map(FILE* in);
		
		// print out matching between places and cells -- note: reading back will not work for files containing more than one map
		//	also for the cell -- place type
		void print_places_map(FILE* out) const;
		
		// copy places_map from a different crecords instance
		void copy_places_map(const crecords& cr) {
			places_maps = cr.places_maps; //assignement operator should be implemented here
			radii = cr.radii;
			pmt = cr.pmt;
		}
		
		// copy places map reversing the direction (i.e. instead of place_id -> cell_id mappings, store cell_id -> place_id mappings)
		void copy_places_map_reverse(const crecords& cr);
		
		//replace cell_id vector according to the given dictionary
		//return: 0: OK, 1: at least one ID was not found
		// note: this only replaces IDs among the records, and not in any of the stored places_maps
		//	it is mainly intended for a spatial aggregation based on a precomputed dictionary of IDs
		int replace_cellids(const hashmap<int,int>& dict);

		//function for better debugging (note: hashmaps do not seem to work well in a debugger)
		static int places_map_count(const hashmap<int,hashset<int> >& places_map, int place_id);
		static int places_map_count(const hashmap<int,hashset<int> >& places_map, int place_id, int cell_id);
		
		void places_map_delete() {
			places_maps.clear();
			radii.clear();
		}
		
		// set only the default places_map explicitely to the radius specified here
		int set_default_places_map(double radius);
		
		
		/******************************************************************
		 * 4. functions for retrieving elements from the CDR records arrays
		 ******************************************************************/
		uint64_t get_size() const { return nrecords; }
		inline uint64_t getnrecords() const { return nrecords; }
		
		//get user ID or return false if i is out of range
		//note: this can be used to loop through all records like this:
		//	for(uint64_t i=0, int64_t uid=0; cr.get_uid(i,&uid); i++) { /* do stuff with uid */ }
		inline bool get_uid(uint64_t i, int64_t* uid) const {
			if(i < uids.size()) { *uid = uids[i]; return true; }
			else return false;
		}
		
		//get timestamp or return false if i is out of range
		//note: this can be used to loop through all records like this:
		//	for(uint64_t i=0, int64_t uid=0; cr.get_uid(i,&uid); i++) { /* do stuff with uid */ }
		inline bool get_ts(uint64_t i, uint32_t* ts1) const {
			if(i < ts.size()) { *ts1 = ts[i]; return true; }
			else return false;
		}
		
		//get record or return false if i is out of range similarly to previous function
		// this could be replaced by a using const iterator
		// this can be used similarly to a C# enumerator interface, e.g.:
		// for(uint64_t i=0,record r; cr.get_record(i,&r); i++) { /* do something with r */ }
		inline bool get_record(uint64_t i, record* r) const {
			if(i < uids.size()) {
				r->uid = uids[i];
				r->ts = ts[i];
				r->cell_id = cells[i];
				if(startstop.size() > 0) {
					if(startstop[i]) r->startstop = 1;
					else r->startstop = 0;
				}
				else r->startstop = -1;
				return true;
			}
			else return false;
		}
		
		// C++ container-like interface
		typedef crecords_iterator iterator;
		typedef crecords_iterator const_iterator;
		typedef ptrdiff_t difference_type;
		typedef size_t size_type;
		typedef record value_type;
		typedef record* pointer;
		typedef record& reference;
		
		const_iterator begin() const;
		const_iterator cbegin() const;
		
		// range-like extension
		typedef crecords_iterator_sentinel sentinel;
		sentinel end() const;
		sentinel cend() const;
		
		// get the first record with the given uid (or nrecords if not found)
		uint64_t get_uid_first(int64_t uid) const;
		
		// get the last+1 record with the given uid (or nrecords if not found -- i.e. that is not necessarily an error, but the return value
		//		of get_uids_first() should be checked for that first)
		uint64_t get_uid_last(int64_t uid) const;
		
		// get an iterator starting at the records of the given user ID
		const_iterator get_uid_iterator_begin(int64_t uid) const;
		
		// get a sentinel for the last record of a given user ID
		sentinel get_uid_iterator_end(int64_t uid) const;
		
		
		/**************************************************************
		 * 5. functions to ensure that data is properly sorted
		 **************************************************************/
		bool use_heapsort; //switch to prefer the new heapsort implementations over the old quicksort ones
		bool use_combined_sort; //use combined quicksort / heapsort with random pivots
		// void set_random_seed_for_sort(uint32_t seed) { rg.init_genrand(seed); }
		void sort_by_uid();
		void sort_by_time();
		void sort_by_uid_time();
		void sort_by_cellid_time();
		void sort_by_uid_cellid();
		void sort_by_cellid_time_uid();
		void sort_by_uid_time_cellid();
		
		sortorder get_sortorder() { return sorted; }
		// check if two crecords classes are sorted by the same -- this is a hack needed by sertest.cpp to compare text and binary data
		//	when explicit sort was requested by the user
		bool is_sorted_same(const crecords& cr2) { return sorted == cr2.sorted; }
		
		void create_cell_idx(); //fill up cell_idx, so that searching by cell_id is made easier
		
		void create_uid_idx(); //fill up uids_idx, so that searching by uid is made easier
		
		
		/*************************************************************
		 * 6. functions calculating statistics about users
		 *************************************************************/
		// note: under normal use, we expect that there are less than 2^32 users and each user has less than 2^32 records
		// 	on the other hand, the total number of records in this class can be more than 2^32
		// to reflect this, the counts are uint32_t types, while the array indices for the records are uint64_t types (or size_t which
		//	is 64-bit on modern compilers to support counts over 2^32); elements in the user index hashmap are also uint64_t, although
		//	the second one could be uint32_t, but it was already used as uint64_t to generate the binary files; this could be changed and
		//	the previous binary files regenerated or converted
		// for now, we implicitly cast to uint32_t here and check for overflow and report error if needed
		
		// return the number of distinct uids
		// note: this recreates the uids_idx hashmap if necessary
		uint64_t get_uids_count();
		
		// get a vector with all the user ids among the records
		// the result must be freed with delete by the caller later
		std::vector<int64_t>* get_uids_p(bool sort = false);
		// same with returning a (local) vector object (using a move constructor in C++11)
		std::vector<int64_t> get_uids(bool sort = false);
		// function doing the actual work for the previous two to fill up a user-supplied vector
		void get_uids(std::vector<int64_t>* r, bool sort = false);
		
		// get a vector with the counts of records for each user
		std::vector<std::pair<int64_t,uint32_t> >* get_uids_cnts_p(bool sort = false);
		// same with returning a (local) vector object (using a move constructor in C++11)
		std::vector<std::pair<int64_t,uint32_t> > get_uids_cnts(bool sort = false);
		// function doing the actual work for the previous two to fill up a user-supplied vector
		void get_uids_cnts(std::vector<std::pair<int64_t,uint32_t> >* r, bool sort = false);
		
		// get the distribution of number of records per user (i.e. there are r[i].second number of users with r[i].first number of matches)
		std::vector<std::pair<uint32_t,uint32_t> >* get_records_dist_p();
		std::vector<std::pair<uint32_t,uint32_t> > get_records_dist();
		void get_records_dist(std::vector<std::pair<uint32_t,uint32_t> >* r);
		
		// cluster the records of a user (gived by uid) based on their coordinates
		// parameters:
		//	tf: optional temporal filter (only include records in the intervals specified there)
		//	res: vector to fill with the results; the contents are the central pair and the number of points in that cluster
		//	cmin: minimum number of points in a cluster to store in the result
		//	rmax: maximum number of clusters to store in the result (zero for these means no limit)
		// return value: 0: OK
		//	> 0: error (wrong sort order, user was not found)
		//	note: if zero points match, or zero cluster are produced that is not reported as an error but results in adding no elements to the result vector
		int cluster_user_coords(int64_t uid, const tfilter* tf, std::vector<std::pair<int,unsigned int> >& res, unsigned int cmin, unsigned int rmax) const;
		
		// find the position and number of records of one user
		// it is an error to call this function with a nonexistent user ID
		std::pair<uint64_t,uint64_t> find_user(int64_t uid) const;
		
//auxilliary functions
	protected:
		
		
		// functions for sorting the array
		
		//swap two records
		inline void pt1_swap(uint64_t i, uint64_t j) {
			int64_t tmp1 = uids[i]; uids[i] = uids[j]; uids[j] = tmp1;
			uint32_t tmp2 = ts[i]; ts[i] = ts[j]; ts[j] = tmp2;
			int tmp3 = cells[i]; cells[i] = cells[j]; cells[j] = tmp3;
		}
		
		// new implementation: only the comparer function has different variants (based on the number of arrays used to compare), defined as overloads
		// there is only one implementation of the sort functions, which (hopefully) call the right version
		template<class T1, class T2, class T3> inline int pt3_cmp(const std::vector<T1>* v1, const std::vector<T2>* v2, const std::vector<T3>* v3,
				uint64_t i, uint64_t j) {
			if( (*v1)[i] < (*v1)[j] ) return 1;
			if( (*v1)[i] > (*v1)[j] ) return -1;
			if( (*v2)[i] < (*v2)[j] ) return 1;
			if( (*v2)[i] > (*v2)[j] ) return -1;
			if( (*v3)[i] < (*v3)[j] ) return 1;
			if( (*v3)[i] > (*v3)[j] ) return -1;
			return 0;
		}
		
		template<class T1, class T2> inline int pt3_cmp(const std::vector<T1>* v1, const std::vector<T2>* v2,
				const std::vector<std::nullptr_t>* v3, uint64_t i, uint64_t j) {
			if( (*v1)[i] < (*v1)[j] ) return 1;
			if( (*v1)[i] > (*v1)[j] ) return -1;
			if( (*v2)[i] < (*v2)[j] ) return 1;
			if( (*v2)[i] > (*v2)[j] ) return -1;
			return 0;
		}
		
		template<class T1> inline int pt3_cmp(const std::vector<T1>* v1, const std::vector<std::nullptr_t>* v2,
				const std::vector<std::nullptr_t>* v3, uint64_t i, uint64_t j) {
			if( (*v1)[i] < (*v1)[j] ) return 1;
			if( (*v1)[i] > (*v1)[j] ) return -1;
			return 0;
		}
		
		template<class T1, class T2 = std::nullptr_t, class T3 = std::nullptr_t>
				void qst3(std::vector<T1>* v1, uint64_t s, uint64_t e, std::vector<T2>* v2 = 0, std::vector<T3>* v3 = 0) {
			if(s >= e) return;
			if(e-s == 1) return;
			if(e-s == 2) {
				if(pt3_cmp(v1,v2,v3,s,s+1) < 0) pt1_swap(s,s+1);
				return;
			}
			
			uint64_t p = s + (e-s)/2;
			uint64_t i;
			pt1_swap(p,e-1);
			p = s;
			for(i=s;i<e-1;i++) {
				int a = pt3_cmp(v1,v2,v3,e-1,i);
				if(a == -1 || (a == 0 && p<i-p) ) {
					if(i>p) pt1_swap(i,p);
					p++;
				}
			}
			if(p < e-1) pt1_swap(p,e-1);
			qst3<T1,T2,T3>(v1,s,p,v2,v3);
			qst3<T1,T2,T3>(v1,p+1,e,v2,v3);
		}
		
		
		// in-place heapsort
		// code originally stolen and adapted from here: http://geeksquiz.com/heap-sort/
		// To heapify a subtree rooted with node i which is
		// an index in arr[]. n is size of heap
		
		// sort by one key
		template <class T1, class T2 = std::nullptr_t, class T3 = std::nullptr_t> void heapify3(std::vector<T1>* v1, uint64_t s, uint64_t n,
				uint64_t i, std::vector<T2>* v2 = 0, std::vector<T3>* v3 = 0)
		{
		    uint64_t largest = i;  // Initialize largest as root
		    uint64_t l = 2*i + 1;  // left = 2*i + 1
		    uint64_t r = 2*i + 2;  // right = 2*i + 2
		 
		    // If left child is larger than root
			if (l < n && pt3_cmp(v1,v2,v3,s+l,s+largest) < 0) largest = l;
			// If right child is larger than largest so far
		    if (r < n && pt3_cmp(v1,v2,v3,s+r,s+largest) < 0) largest = r;
		 
		    // If largest is not root
		    if (largest != i)
		    {
		        pt1_swap(s+i,s+largest); //make the root the largest
		        // Recursively heapify the affected sub-tree
		        heapify3(v1, s, n, largest, v2, v3);
		    }
		}
		 
		// main function to do heap sort
		template <class T1, class T2 = std::nullptr_t, class T3 = std::nullptr_t> void heapSort3(std::vector<T1>* v1, uint64_t s, uint64_t e,
				std::vector<T2>* v2 = 0, std::vector<T3>* v3 = 0)
		{
			uint64_t n = e-s;
			if(n <= 1) return;
		    // Build heap (rearrange array)
		    uint64_t i = n/2-1;
		    while(1) { heapify3<T1,T2,T3>(v1, s, n, i, v2, v3); if(i==0) break; i--; }
		    
		    // One by one extract an element from heap
		    i = n-1;
		    for (i=n-1; i>0; i--)
		    {
		        // Move current root to end
		        pt1_swap(s, s+i);
		        // call max heapify on the reduced heap
		        heapify3<T1,T2,T3>(v1, s, i, 0, v2, v3);
		    }
		}
		
		template<class T1, class T2 = std::nullptr_t, class T3 = std::nullptr_t>
				void qst3c(std::vector<T1>* v1, uint64_t s, uint64_t e, std::vector<T2>* v2 = 0, std::vector<T3>* v3 = 0, unsigned int depth = 0) {
			if(s >= e) return;
			if(e-s == 1) return;
			if(e-s == 2) {
				if(pt3_cmp(v1,v2,v3,s,s+1) < 0) pt1_swap(s,s+1);
				return;
			}
			
			if(e-s < 12288 || depth > 64) { //manual threshold, for these, use the heapsort which has a better worst-case performance and is only a little bit slower
					// on average for this size -- 12288 means 192KiB, which should fit in the L2 cache in most modern CPUs
					// + extra threshold: depth is the recursion depth, 64 is already huge (i.e. log_2(500M) ~ 29) to protect from
					//	stack overflow (on debug builds where the recursion is not eliminated) and O(N^2) behavior which happens too often
				heapSort3<T1,T2,T3>(v1,s,e,v2,v3);
				return;
			}
			
			uint64_t p = s + (e-s)/2;
			uint64_t i;
			pt1_swap(p,e-1);
			p = s;
			for(i=s;i<e-1;i++) {
				int a = pt3_cmp(v1,v2,v3,e-1,i);
				if(a == -1 || (a == 0 && p<i-p) ) {
					if(i>p) pt1_swap(i,p);
					p++;
				}
			}
			if(p < e-1) pt1_swap(p,e-1);
			qst3c<T1,T2,T3>(v1,s,p,v2,v3,depth+1);
			qst3c<T1,T2,T3>(v1,p+1,e,v2,v3,depth+1);
		}
	 
		
		// binary search in the ts array (note: array needs to be sorted by time, this is currently ensured by the read functions
		//		and the functions calling this one)
		// return value:
		//		if start == true, the first occurrence of t1 or the smallest ts value > t1 if it is not found
		//			or e if t1 > all ts values
		//		if start == false, the last occurrence of t1 or the largest ts value < t1 if it is not found
		//			or e if t1 < all ts values (note: 0 would mean that ts[0] <= t1)
		// this was a range can be found with searching for tsearch(start,true) and tsearch(end,false);
		// limit search between s and e (s inclusive, e exclusive; if omitted, search the whole range)
		uint64_t tsearch(uint32_t t1, bool start, uint64_t s, uint64_t e) const;
		
		uint64_t tsearch(uint32_t t1, bool start) const {
			return tsearch(t1,start,0,nrecords);
		}
		
		
		// fill the uids_idx hashmap with the number occurrences of each user
		void uids_count();
		
		// compare record i in this crecords class with records other_i ... e2 in the other class
		//	(helper function for the compare_users_onematch() above)
		// return value:
		//	1: other_i is a match
		//	0: other_i is not a temporal match
		//	-1: there is an impossible match (temporal match but not spatial) in other_i ... e2
		//	-2: error (place_id not found)
		//~ int compare_users_onematch_r(uint64_t i, unsigned int dt, const crecords& other_records, uint64_t other_i,
				//~ uint64_t e2, bool ignore_missing, bool use_other_map) const;
		

#ifdef USECEREAL
	private:
		
		friend class cereal::access;
		
		// note: load should only be called with an empty class
		template <class Archive> void load(Archive& ar, const uint32_t version) {
			if(version >= 4) // the Voronoi-tessellations, place coordinates and place-cell mappings are not saved anymore
				ar(
						uids,
						ts,
						cells,
						nrecords,
						startstop,
						
						sorted,
						uids_idx,
						cells_idx,
						uids_idx_cnts
					);
			
			else throw new std::runtime_error("crecords::load(): unsupported binary file format!\n");
		}
		
		template <class Archive> void save(Archive& ar, const uint32_t version) const {
			ar(uids,ts,cells,nrecords,startstop,sorted,uids_idx,cells_idx,uids_idx_cnts);
		}
#endif

	public:
		//(for testing purposes)
		//compare the data in this class to an other one (implemented in sertest.cpp)
		//	note: if ignore_places == true, the places_map and all related fields are ignored
		bool compare_crecords(crecords& other, bool ignore_places = false);
		
		//read data from a specified file
		static int read_crecords_serialized(crecords& cr, const char* fn, bool zip = false);
		
		//write data to a specified file (note: compression is not supported here yet, could be easily done with gzip after writing)
		static int write_crecords_serialized(const crecords& cr, const char* fn);
};


#ifdef USECEREAL
CEREAL_CLASS_VERSION(crecords,4);
#endif


// helper function for generating a random number in a given range using std::random
// R should be a random number generator returning T
// the value returned is in the range [0,max-1]
// result is undefined if max == 0
template<class R> typename R::result_type get_random_range(R& rg, typename R::result_type max) {
	std::uniform_int_distribution<typename R::result_type> dist(0,max-1);
	return dist(rg);
	/*
	
	typename R::result_type res;
	typename R::result_type rmin = rg.min();
	typename R::result_type rmax = rg.max() - rmin;
	typename R::result_type rem = rmax % max + 1;
	if(rem == max) rem = 0;
	typename R::result_type ign = rmax - rem;
	do {
		res = rg() - rmin;
	} while(res > ign);
	return res % max;*/
}


/**********************************************************
 * separate class for storing the temporal distribution
 * of points from one of the datasets
 * can be generated from the crecords class
 **********************************************************/
class tdist {
	private:
		//minimum and maximum time (as UNIX timestamp, inclusive)
		uint32_t tmin;
		uint32_t tmax;
		
		//array holding the distribution with 1 second precision; size: tmax - tmin + 1,
		//	sdist[i] corresponds to the number of records at timestamp tmin+i
		//	(note: for a one week dataset, this means 604800 elements, i.e. ~2.3MB)
		std::deque<uint32_t> sdist;
		
		//array holding the distribution with 1 minute precision (to speed up computations), size: (tmax - tmin) / 60 + 1
		//	mdist[i] corresponds to the number of records between timestamps tmin + i*60, tmin + i*60 + 59 (inclusive)
		//	(a records with timestamp t goes into the bin of (t-tmin)/60
		std::deque<uint32_t> mdist;
		
		// cumulative distributions (used for generating a random timestamp based on the empirical distributions)
		// these are only generated when needed
		std::vector<uint64_t> scdf; // note: we need to use std::vector to be able to use std::lower_bound() for binary search
		std::mt19937_64 rg; // random number generator
		
		uint64_t nrecords; //total number of records stored here
	
	public:
		tdist() {
			tmin = 0;
			tmax = 0;
			nrecords = 0;
		}
		
		// create temporal distribution based on the supplied crecords class
		// note: if tsorted == true, the records are assumed to be sorted by time (so that we do not separately go over all values to find the
		//		minimum and maximum)
		// return values: 0: OK; 1: error allocating memory (note: depending on the compiler and library, this might throw an exception in that case too)
		template<class iterator, class sentinel>
		void createdist(iterator it, const sentinel& end) {
			// clear any previous data
			mdist.clear();
			sdist.clear();
			tmin = std::numeric_limits<int>::max(); //note: this is 2^32-1
			tmax = 0;
			nrecords = 0;
			uint32_t ts;
			
			// determine and add first element
			if(it == end) return;
			ts = (*it).ts;
			tmin = ts;
			tmax = ts;
			uint32_t tmin1 = tmin - tmin%604800; // position one week along the given timestamp
			uint32_t tmax1 = tmin1 + 604799;
			
			// allocate arrays
			sdist.resize(tmax1-tmin1+1,0); // this will be always 604800, i.e. exactly one week
			mdist.resize( (tmax1-tmin1)/60 + 1 ); // this will be always 10080
			
			// go through all records, classify them into the bins
			// do not use the supplied iterator, so iteration can be restarted if not sorted
			for(;it!=end;++it) {
				ts = (*it).ts;
				if(ts < tmin) tmin = ts;
				if(ts > tmax) tmax = ts;
				
				if(ts < tmin1) { // resize the arrays -> add tmin1-ts elements (round up to one day increments)
					uint32_t tmin2 = ts - ts%86400;
					sdist.insert(sdist.begin(),tmin1-tmin2,0);
					mdist.insert(mdist.begin(),(tmin1-tmin2)/60,0); // note: (tmin1-tmin2)%60 == 0 always, since both are day boundaries
					tmin1 = tmin2;
				}
				if(ts > tmax1) { // add ts-tmax1 elements at the end (round up to one day increments again)
					uint64_t tmax2 = ts + (86400 - ts%86400); // note: tmax2%86400 == 0 now, this is one more than the real end
					sdist.insert(sdist.end(),tmax2-tmax1-1,0); // -1 to make it an exact multiple of 86400
					mdist.insert(mdist.end(),(tmax2-tmax1-1)/60,0); // (tmax2-tmax1-1)%60 == 0 in all cases
					tmax1 = tmax2-1;
				}
				
				sdist[ts-tmin1]++;
				mdist[(ts-tmin1)/60]++;
				nrecords++;
			}
			
			// resize the arrays back, so that tmin1 == tmin and tmax1 == tmax
			if(tmin > tmin1) {
				// delete the first elements
				sdist.erase(sdist.begin(), sdist.begin() + (tmin - tmin1) );
				if(tmin - tmin1 >= 60) {
					// elements from mdist need to be deleted as well
					uint32_t tmin2 = tmin1 - tmin1%60;
					mdist.erase(mdist.begin(), mdist.begin() + (tmin2-tmin1)/60);
				}
			}
			if(tmax1 > tmax) {
				// delete the last elements
				sdist.erase(sdist.begin() + (tmax - tmin + 1), sdist.end());
				if(tmax1 >= tmax + 60) {
					// delete from mdist as well
					uint32_t mlast = tmax/60 - tmin/60;
					mdist.erase(mdist.begin() + mlast + 1, mdist.end());
				}
			}
		}

		
		inline uint32_t get_tmax() const { return tmax; }
		inline uint32_t get_tmin() const { return tmin; }
		
		inline uint64_t getnrecords() const { return nrecords; } // total number of records
		inline uint64_t getnt(uint32_t t1, uint32_t t2) const { // get the number of records between t1 and t2 (inclusive)
			if(t2 > tmax) t2 = tmax;
			if(t1 < tmin) t1 = 0;
			else t1 -= tmin;
			if(t2 < tmin) return 0;
			else t2 -= tmin;
			uint32_t t = t1;
			uint64_t r = 0;
			for(;t%60 && t<=t2;t++) r += sdist[t];
			uint32_t m2 = t2/60;
			uint32_t m = t/60;
			for(;m<m2;m++) r += mdist[m];
			t = m*60;
			for(;t<=t2;t++) r +=  sdist[t];
			return r;
		}
		inline double getpt(uint32_t t1, uint32_t t2) const { // get the probability that a record is between t1 and t2 (inclusive)
			// this is just the number of records there / all records, so this is mainly a convenienve function
			return ((double)getnt(t1,t2)) / ((double)nrecords);
		}
		
		// get the number of points in the ranges [t-dt,t+dt]
		// for nonoverlapping ranges, calling getnt() for each of the ranges might be more efficient; here extra care is taken to count records
		//	from the overlapping parts only once
		// use the first s elements in the given vector only
		// the sorted parameter is a hint whether the t vector is already sorted (if not, it is sorted by std::sort();
		//	note that this modifies the original vector!)
		// note: additions and subtractions (for t[i]-dt, t[i]+dt) are assumed not to overflow, otherwise the result is undefined and possibly a segfaults
		// note: instead of a vector, a more general iterator could be used for sorted inputs, i.e. an iterator returning elements from a crecords class
		inline uint64_t getnt2(std::vector<uint32_t>* t0, size_t s, uint32_t dt, bool sorted = false) const {
			if(!sorted) std::sort(t0->begin(),t0->begin() + s);
			uint64_t r = 0;
			const std::vector<uint32_t>& t = *t0; // for convenience
			if(s == 0) return 0;
			if(s == 1) return getnt(t[0]-dt,t[0]+dt);
			uint32_t t1 = t[0]-dt;
			uint32_t t2 = t[0]+dt;
			size_t i = 1;
			while(1) {
				for(;i < s;i++) { //merge overlapping intervals
					if(t[i]-dt > t2+1) break;
					t2 = t[i]+dt;
					//assure that intervals are really sorted
					if(t[i]-dt < t1) return getnt2(t0,dt,false); //not really sorted, re-sort and re-try
				}
				r += getnt(t1,t2);
				if(i == s) break;
				t1 = t[i]-dt;
				t2 = t[i]+dt;
				i++;
			}
			return r;
		}
		
		inline uint64_t getnt2(std::vector<uint32_t>* t0, uint32_t dt, bool sorted = false) const {
			return getnt2(t0, t0->size(), dt, sorted);
		}
		
		inline double getpt2(std::vector<uint32_t>* t, size_t s, uint32_t dt, bool sorted = false) const { // same as the previous one, but return probabilities
			return ((double)getnt2(t,s,dt,sorted)) / ((double)nrecords);
		}
		inline double getpt2(std::vector<uint32_t>* t, uint32_t dt, bool sorted = false) const { // same as the previous one, but return probabilities
			return ((double)getnt2(t,t->size(),dt,sorted)) / ((double)nrecords);
		}
		
		// similar to the previous one, but return the number of points not in the intervals (i.e. subtract the result from nrecords)
		inline uint64_t getnt2i(std::vector<uint32_t>* t, size_t s, uint32_t dt, bool sorted = false) const {
			return nrecords - getnt2(t,s,dt,sorted);
		}
		inline uint64_t getnt2i(std::vector<uint32_t>* t, uint32_t dt, bool sorted = false) const {
			return nrecords - getnt2(t,t->size(),dt,sorted);
		}
		
		inline double getpt2i(std::vector<uint32_t>* t, size_t s, uint32_t dt, bool sorted = false) const { // same as the previous one, but return probabilities
			return ((double)getnt2i(t,s,dt,sorted)) / ((double)nrecords);
		}
		inline double getpt2i(std::vector<uint32_t>* t, uint32_t dt, bool sorted = false) const { // same as the previous one, but return probabilities
			return ((double)getnt2i(t,t->size(),dt,sorted)) / ((double)nrecords);
		}
		
		
		// generate random timestamps based on the distributions
		// the given seed is used to initialize the random number generator if useseed == true
		//	(otherwise, the default initialization in the implementation is used)
		// return 0 on success, 1 on error (memory allocation failure or if the PDFs were not initialized previously)
		int rg_prepare(bool useseed = false, uint64_t seed = 0);
		
		// get a random timestamp according to the distribution from the given range
		// if min == max == 0, the whole time range in this instance is used
		// if the range is invalid or the cdfs are uninitialized and exception is thrown
		// if no timestamps are in the range, a timestamp out of the range is returned (0, if it is not in the range)
		// same function, but use the supplied rng, so that it can be made thread-safe
		uint32_t get_random_ts_r(std::mt19937_64& rg1, uint32_t min = 0, uint32_t max = 0) const;
		inline uint32_t get_random_ts(uint32_t min = 0, uint32_t max = 0) {
			return get_random_ts_r(rg,min,max);
		}
		
		
};


// "iterator" class generating random timestamps using the previous one for filling up a crecords with records with random timestamps
// note that this only generates timestamps, the user ID and the cell ID in the returned struct is always set to 0
class tdistriterator : std::iterator<std::input_iterator_tag, record> {
	private:
		uint64_t i; // counter (for iterator comparison)
		tdist* td; // distribution to generate timestamps from
		uint32_t tmin; // parameters for generating random timestamps
		uint32_t tmax;
		uint32_t tsnext; // store the next timestamp, so operator*() always returns the same, while operator++() generates the next timestamp
		// note: this means that the last generated value will generally be wasted (if used in a usual iterator-loop)
		tdistriterator() { td = 0; i = 0; rg = 0; }
		std::mt19937_64* rg;
		
		void generate_next() {
			if(td == 0) throw new std::runtime_error("tdistriterator(): invalid iterator used!\n");
			if(rg) tsnext = td->get_random_ts_r(*rg,tmin,tmax);
			else tsnext = td->get_random_ts(tmin,tmax);
			if( (tmin == 0 && tmax == 0 && tsnext == 0) || ( (tmin>0 || tmax>0) && (tsnext < tmin || tsnext > tmax) ) ) {
				throw new std::runtime_error("tdistriterator::generate_next(): error generating random timestamps!\n");
			}
			
		}
	public:
		tdistriterator(tdist* t, uint32_t tmin_ = 0, uint32_t tmax_ = 0) { i = 0; td = t; tmin = tmin_; tmax = tmax_; generate_next(); }
		tdistriterator(tdist* t, std::mt19937_64* rg_, uint32_t tmin_ = 0, uint32_t tmax_ = 0) { i = 0; td = t; tmin = tmin_; tmax = tmax_; rg = rg_; generate_next(); }
		tdistriterator(uint64_t N) { i = N; td = 0; rg = 0; tsnext = 0; } // create an invalid iterator, only useful for comparing to the end
		record operator *() const {
			record r;
			r.uid = 0;
			r.ts = tsnext;
			r.cell_id = 0;
			r.startstop = -1;
			return r;
		}
		void operator++() { generate_next(); i++; }
		bool operator==(const tdistriterator& r) const { return i == r.i; }
		bool operator!=(const tdistriterator& r) const { return i != r.i; }
		friend class tdistriterator_sentinel;
		friend class tdistrpairiterator;
};


// event length distribution (used by the next class)
class tlendistr {
	private:
		std::vector<uint32_t> lengths; // possible lengths
		std::vector<uint64_t> lcdf; // CDF corresponding to the previous (including raw frequencies, so that 64-bit numbers are generated for sampling
	public:
		tlendistr() {  }
		// initialize the length distributions (read from the given file)
		int init_length_distr(FILE* f);
		friend class tdistrpairiterator;
};

// extension of the previous iterator so that records are sometimes generated in pairs according to a separate length distribution
//	(read from a file separately)
class tdistrpairiterator : std::iterator<std::input_iterator_tag, record> {
	private:
		uint64_t i; // counter (for comparison)
		tdistriterator tdit; // iterator to draw the event timestamps from
		uint32_t tsnext1; // next timestamp to return
		uint32_t tsnext2; // second next timestamp (end of call / trip)
		//~ std::vector<uint32_t> lengths; // possible lengths
		//~ std::vector<uint64_t> lcdf; // CDF corresponding to the previous (including raw frequencies, so that 64-bit numbers are generated for sampling
		const tlendistr* lendistr;
		tdistrpairiterator():i(0),tsnext1(0),tsnext2(0),lendistr(0) {  }
		
		void generate_next() {
			if(tsnext2 > 0) { // call / trip end is stored
				tsnext1 = tsnext2;
				tsnext2 = 0;
			}
			else { // generate a new call / trip based on the length distribution
				tdit.generate_next();
				if( lendistr == 0 || lendistr->lengths.size() == 0 || lendistr->lengths.size() != lendistr->lcdf.size() )
					throw new std::runtime_error("tdistrpairiterator::generate_next(): length distribution was not properly initialized!\n");
				uint64_t max = lendistr->lcdf[lendistr->lcdf.size()-1];
				uint64_t r1 = get_random_range(*(tdit.rg),max);
				size_t r2 = std::lower_bound(lendistr->lcdf.begin(), lendistr->lcdf.end(), r1) - lendistr->lcdf.begin();
				tsnext1 = tdit.tsnext;
				tsnext2 = tsnext1 + lendistr->lengths[r2];
			}
		}
		
		void check_init(const tlendistr* lendistr_) {
			if(lendistr_ == 0) throw new std::runtime_error("tdistrpairiterator: length distribution needs to be given at initialization!\n");
			if( lendistr_->lengths.size() == 0 || lendistr_->lengths.size() != lendistr_->lcdf.size() ) throw new std::runtime_error("tdistrpairiterator: length distribution needs to be initialized!\n");
			if(tdit.rg == 0) throw new std::runtime_error("tdistrpairiterator: random number generator needs to be given at initialization!\n");
			lendistr = lendistr_;
			generate_next();
		}
		
	public:
		tdistrpairiterator(tdistriterator& tdit_, const tlendistr* lendistr_): tdit(tdit_),i(0),tsnext1(0),tsnext2(0) { check_init(lendistr_); }
		tdistrpairiterator(tdistriterator&& tdit_, const tlendistr* lendistr_): tdit(tdit_),i(0),tsnext1(0),tsnext2(0) { check_init(lendistr_); }
		record operator *() const {
			if(tsnext1 == 0) throw new std::runtime_error("tdistrpairiterator: length distribution was not properly initialized!\n");
			record r;
			r.uid = 0;
			r.ts = tsnext1;
			r.cell_id = 0;
			r.startstop = -1;
			return r;
		}
		void operator++() { generate_next(); i++; }
		// note: comparison is deprecated, use the tdistriterator_sentinel class
		friend class tdistriterator_sentinel;
};



// separate sentinel class for these (essentially just the number of records to generate)
class tdistriterator_sentinel {
	private:
		uint64_t i;
	public:
		tdistriterator_sentinel():i(0) {  }
		tdistriterator_sentinel(uint64_t i_):i(i_) {  }
		bool operator==(const tdistriterator& r) const { return i <= r.i; } // signal the end if the other iterator is already above the number here
		bool operator!=(const tdistriterator& r) const { return i > r.i; } 
		bool operator==(const tdistrpairiterator& r) const { return i <= r.i; } // signal the end if the other iterator is already above the number here
		bool operator!=(const tdistrpairiterator& r) const { return i > r.i; } 
};

static inline bool operator==(const tdistriterator& i, const tdistriterator_sentinel& s) { return s.operator==(i); }
static inline bool operator!=(const tdistriterator& i, const tdistriterator_sentinel& s) { return s.operator!=(i); }
static inline bool operator==(const tdistrpairiterator& i, const tdistriterator_sentinel& s) { return s.operator==(i); }
static inline bool operator!=(const tdistrpairiterator& i, const tdistriterator_sentinel& s) { return s.operator!=(i); }

const static tdistriterator tdistriterator_null( (uint64_t) 0 );

// example usage
static inline void crecords_fill_random_ts(crecords* cr, tdist* td, uint64_t N, uint32_t tmin = 0, uint32_t tmax = 0) {
	tdistriterator r(td,tmin,tmax);
	cr->load_records_custom(r,tdistriterator_sentinel(),N,true);
}



// "iterator" class to fill up crecords with totally random data (generated using lrand48() from the C standard library
// this is not thread-safe in it's current form !!
class randomiterator : std::iterator<std::input_iterator_tag, record> {
	private:
		uint64_t i;
	public:
		randomiterator() { i = 0; }
		randomiterator(uint64_t N) { i = N; }
		record operator *() {
			record r;
			r.uid = lrand48();
			r.ts = (uint32_t)lrand48();
			r.cell_id = (int)lrand48();
			r.startstop = -1;
			return r;
		}
		void operator++() { i++; }
		bool operator==(const randomiterator& r) const { return i == r.i; }
		bool operator!=(const randomiterator& r) const { return i != r.i; }
};

// example usage
static inline void crecords_fill_random_data(crecords* cr, uint64_t N) {
	randomiterator r;
	cr->load_records_custom(r,randomiterator(N),N);
}



// helper class to store the number of records per user and pick a user randomly in a range
// note: this requires C++11
class user_rand {
	protected:
		std::vector<std::pair<int64_t,uint32_t> > uids; //user ids and counts, sorted by counts (second)
		uint32_t nusers; // number of users (<2^32)
		uint32_t cntmax; // maximum count in uids
		uint32_t cntmin; // minimum count in uids
		void sort_uids();
		std::mt19937 r;
		
	public:
		user_rand() { nusers = 0; cntmax = 0; cntmin = 0; }
		user_rand(std::vector<std::pair<int64_t,uint32_t> > uids1) { uids = uids1; sort_uids(); }
		user_rand(std::vector<std::pair<int64_t,uint32_t> > uids1, uint32_t seed) { uids = uids1; sort_uids(); r.seed(seed); }
		void set_uids(std::vector<std::pair<int64_t,uint32_t> > uids1) { uids = uids1; sort_uids(); }
		
		inline void init_genrand(uint32_t s) { r.seed(s); }
	//	inline void init_genrand_array(uint32_t* init_key, int key_length) { r.init_by_array(init_key, key_length); }
		
		uint32_t getnusers() { return nusers; }
		uint32_t getmin() { return cntmin; }
		uint32_t getmax() { return cntmax; }
		
		// return a random user id with record count between min and max (inclusive), or -1 if no such user exists
		// if min == max == 0, just pick any of the users at random
		// if nrecords is given, return the number of records the chosen user actually has in it
		int64_t get_random_uid(uint32_t min, uint32_t max, uint32_t* nrecords = 0);
		inline int64_t get_random_uid_all(uint32_t* nrecords = 0) { return get_random_uid(0,0,nrecords); }
};



// c++ wrapper for reading from FILE* streams (e.g. the one returned by the above functions)
#include <iostream>
#include <streambuf>
class stdiobuf
    : public std::streambuf
{
private:
    FILE* d_file;
    char  d_buffer[8192];
public:
    stdiobuf(FILE* file): d_file(file) {}
    int underflow() {
        if (gptr() == egptr() && d_file) {
            size_t size = fread(d_buffer, 1, 8192, d_file);
            setg(d_buffer, d_buffer, d_buffer + size);
        }
        if(gptr() == egptr()) return traits_type::eof();
        else return traits_type::to_int_type( * (gptr()) );
    }
};

class stdiostream {
	private:
		stdiobuf buf;
		std::istream is;
	public:
		stdiostream(FILE* f): buf(f), is(&buf) {  }
		inline std::istream& stream() { return is; }
};


/*
 * C++ iterators for iterating over a crecords class or reading data from a text file
 * 	+ iterator which filters based on user IDs, cell IDs or timestamps
 * these can be chained and input to add_records_custom
 */
class crecords_iterator : std::iterator<std::input_iterator_tag, record> {
	protected:
		bool is_end; // flag to indicate if we are already at the end
		const crecords* cr; // records to iterate from
		uint64_t i; // counter (position in crecords)
		crecords_iterator() { is_end = true; cr = 0; i = 0; }
		record r;
		
	public:
		crecords_iterator(const crecords* cr_) {
			cr = cr_;
			is_end = false;
			i = 0;
			if(cr->getnrecords() == 0) is_end = true;
			else if(!cr->get_record(0,&r)) is_end = true;
		}
		crecords_iterator(const crecords* cr_, uint64_t i_) {
			cr = cr_;
			is_end = false;
			i = i_;
			if(i >= cr->getnrecords()) is_end = true;
			else if(!cr->get_record(i,&r)) is_end = true;
		}
		
		record operator *() const {
			if(is_end) throw new std::runtime_error("crecords_iterator(): iterator used after reaching the end!\n");
			return r;
		}
		const record* operator ->() const {
			if(is_end) throw new std::runtime_error("crecords_iterator(): iterator used after reaching the end!\n");
			return &r;
		}
		void operator++() {
			i++;
			if(!cr->get_record(i,&r)) is_end = true;
		}
		bool operator==(const crecords_iterator& r) const { return is_end == r.is_end; }
		bool operator!=(const crecords_iterator& r) const { return is_end != r.is_end; }
		
		static const crecords_iterator it_end() {
			const crecords_iterator end;
			return end;
		}
		
		bool is_it_end() const {
			return is_end;
		}
		friend class crecords_iterator_sentinel;
};

class crecords_iterator_sentinel { // sentinel class for a crecords iterator (compare the index)
	private:
		uint64_t i;
	public:
		crecords_iterator_sentinel():i(0) {  }
		crecords_iterator_sentinel(uint64_t i_):i(i_) {  }
		
		bool operator==(const crecords_iterator& it) const { return it.i >= i || it.is_it_end(); }
		bool operator!=(const crecords_iterator& it) const { return it.i < i && !it.is_it_end(); }
};

static inline bool operator==(const crecords_iterator& it, const crecords_iterator_sentinel& s) { return s.operator==(it); }
static inline bool operator!=(const crecords_iterator& it, const crecords_iterator_sentinel& s) { return s.operator!=(it); }


class tsv_iterator : std::iterator<std::input_iterator_tag, record> {
protected:
	FILE* in;
	uint64_t lines_max;
	uint64_t header_skip;
	uint64_t line;
	bool is_end;
	record r;
	bool read_startstop;
	
	tsv_iterator() { is_end = true; in = 0; line = 0; lines_max = 0; read_startstop = false; }
	
	void read_next() {
		line++;
		if(lines_max && line > lines_max + header_skip) { is_end = true; return; }
		//!! note: %lu might not be portable !!
		//note: we use %*[ \t] instead of white space so that a newline cannot match there -- problem: numeric conversion specifiers skip newlines?
		int a = 0;
		int rs = 0;
		if(read_startstop) {
			int s;
			a = fscanf(in,"%ld%*[ \t]%u%*[ \t]%d%*[ \t]%d",&(r.uid),&(r.ts),&(r.cell_id),&s);
			if(s) r.startstop = 1;
			else r.startstop = 0;
			rs = 4;
		}
		else {
			a = fscanf(in,"%ld%*[ \t]%u%*[ \t]%d",&(r.uid),&(r.ts),&(r.cell_id));
			r.startstop = -1;
			rs = 3;
		}
		if(a == EOF) { is_end = true; return; }
		if(a != rs) {
			fprintf(stderr,"tsv_iterator::read_next(): Invalid data in line %lu!\n",line);
			is_end = true;
			return;
		}
		
		do {
			a = getc(in);
		} while( ! (a == '\n' || a == EOF) );
	}
	
public:
	tsv_iterator(FILE* in_, uint64_t header_skip_ = 0, uint64_t lines_max_ = 0, bool read_startstop_ = false) {
		in = in_;
		header_skip = header_skip_;
		read_startstop = read_startstop_;
		for(uint64_t j=0;j<header_skip;j++) {
			int a;
			do {
				a = getc(in);
			} while( ! (a == '\n' || a == EOF) );
			if(a == EOF) break;
		}
		if(feof(in) || ferror(in)) is_end = true;
		line = header_skip;
		lines_max = lines_max_;
		read_startstop = read_startstop_;
		is_end = false;
		read_next();
	}
	
	record operator *() {
		if(is_end) throw new std::runtime_error("tsv_iterator(): iterator used after reaching the end!\n");
		return r;
	}
	void operator++() {
		read_next();
	}
	bool operator==(const tsv_iterator& r) const { return is_end == r.is_end; }
	bool operator!=(const tsv_iterator& r) const { return is_end != r.is_end; }
		
	static const tsv_iterator it_end() {
		const tsv_iterator end;
		return end;
	}
	
	bool is_it_end() const {
		return is_end;
	}
};


// iterator to filter results (obtained via an iterator returning records) based on user ID
template <class it, class sentinel>
class filter_iterator : std::iterator<std::input_iterator_tag, record> {
	private:
		bool is_end; // flag to indicate if we are already at the end
		it current;
		const sentinel end;
		record r; // current record
		const std::unordered_set<int64_t>* users;
		const std::unordered_set<int>* cells;
		uint32_t tmin;
		uint32_t tmax;
		
		filter_iterator(); // no default constructor, always has to be given a range to filter
			// for "sentinels" (filter_iterators denoting the end of a range) use either record_iterator_sentinel
			// or the get_end function of a previously constructed filter_iterator
		
		bool check_r() {
			bool found = true;
			if(tmax > 0 && tmax > tmin) if(r.ts < tmin || r.ts >= tmax) found = false;
			if(found && users && users->count(r.uid) == 0) found = false;
			if(found && cells && cells->count(r.cell_id) == 0) found = false;
			return found;
		}
		
		void it_find_next() {
			for(is_end=true;current != end && is_end == true;++current) {
				r = *current;
				bool found = check_r();
				if(found) is_end = false; // note: loop will break after advancing current
			}
		}
		
		
	public:
		filter_iterator(it current_, const sentinel end_, const std::unordered_set<int64_t>* users_ = 0, const std::unordered_set<int>* cells_ = 0,
				uint32_t tmin_ = 0, uint32_t tmax_ = 0):current(current_),end(end_) {
			users = users_;
			cells = cells_;
			if(tmin_ > 0 || tmax_ > 0) set_tsfilter(tmin,tmax);
			is_end = false;
			it_find_next(); // find the first matching user -- contrary to C#, where an iterator has a status indicating whether it reached the end
				// (and "dereferencing" is actually a member function returning bool indicating if it has already reached the end),
				// in C++, the dereference operator (operator *()) must always return a valid result and checking for the end is done with comparing
				// to an other instance of the iterator -> this means that we have to check if there will be any results in advance, and calculate
				// the next result at each dereference
		}
		
		void set_tsfilter(uint32_t tmin_, uint32_t tmax_) {
			if(tmax_ <= tmin_) {
				fprintf(stderr,"filter_iterator::set_tsfilter(): invalid time interval (%u,%u)!\n",tmin_,tmax_);
				throw new std::runtime_error("filter_iterator::set_tsfilter(): invalid time interval");
			}
			tmin = tmin_;
			tmax = tmax_;
			if(!is_end && check_r() == false) it_find_next();
		}
		void set_userfilter(const std::unordered_set<int64_t>* users_) {
			users = users_;
			if(!is_end && check_r() == false) it_find_next();
		}
		void set_cellsfilter(const std::unordered_set<int>* cells_) {
			cells = cells_;
			if(!is_end && check_r() == false) it_find_next();
		}
		
		record operator *() {
			if(is_end) throw new std::runtime_error("filter_iterator(): iterator used after reaching the end!\n");
			return r;
		}
		void operator++() {
			it_find_next(); // advance i to the next possible position (or set is_end)
		}
		bool operator==(const filter_iterator& r) const { return is_end == r.is_end; }
		bool operator!=(const filter_iterator& r) const { return is_end != r.is_end; }
		
		const filter_iterator get_end() {
			const filter_iterator end_(current,end);
			end_.is_end = true;
			return end_;
		}
		
		bool is_it_end() const {
			return is_end;
		}
};


// sentinel class for any of the previous three iterators (crecords, tsv and filter)
class record_iterator_sentinel {
	// empty class just used in the comparison functions
};

// "comparisons" which will work for any of the previous and signal the end
template<class recorditerator>
static inline bool operator==(const record_iterator_sentinel& s, const recorditerator& i) {
	return i.is_it_end(); // i == s if i is at the end (can be used in a C++ style loop)
}
template<class recorditerator>
static inline bool operator==(const recorditerator& i, const record_iterator_sentinel& s) {
	return i.is_it_end(); // i == s if i is at the end (can be used in a C++ style loop)
}

template<class recorditerator>
static inline bool operator!=(const record_iterator_sentinel& s, const recorditerator& i) {
	return !(i.is_it_end()); // i == s if i is at the end (can be used in a C++ style loop)
}
template<class recorditerator>
static inline bool operator!=(const recorditerator& i, const record_iterator_sentinel& s) {
	return !(i.is_it_end()); // i == s if i is at the end (can be used in a C++ style loop)
}

// empty sentinel class (always returns false) to be used with "infinite" iterators, where a different way of stopping is available
class null_sentinel {
	// empty class just used in the comparison functions
};

template<class iterator>
static inline bool operator==(const null_sentinel& s, const iterator& i) {
	return false;
}
template<class iterator>
static inline bool operator==(const iterator& i, const null_sentinel& s) {
	return false;
}

template<class iterator>
static inline bool operator!=(const null_sentinel& s, const iterator& i) {
	return true;
}
template<class iterator>
static inline bool operator!=(const iterator& i, const null_sentinel& s) {
	return true;
}


#endif // CDRRECORDS_H

