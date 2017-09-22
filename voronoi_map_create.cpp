/*
 * voronoi_map_create.cpp
 * creating map between places and cells using pregenerated Voronoi tessellation of cells
 * write out the result to a file for using later with the matching_* programs
 * 
 * note: this is a standalone program, it does not depend on the crecords class or any other functionality in those files
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
 * ./vtest -c ../matching/voronoi/cell_coords.csv -p ../lta/stops.dat -vc ../matching/voronoi/vertices.txt -vv ../matching/voronoi/vertex_to_cell.txt -ve ../matching/voronoi/ridge_to_vertex2.txt -vu ../matching/voronoi/unique_to_cell.txt > ../matching/voronoi/vt500m.out
 * note: runtime ~22s (debug version)
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility> //std::pair

// note: currently, this program requires C++11, but this is not necessary:
//	the C++98 std::map and std::set could be used by replacing the following two lines with appropriate #defines
template <class T, class U> using hashmap = std::unordered_map<T,U>;
template <class T> using hashset = std::unordered_set<T>;



// helper struct to hold all the data 
struct mapper {
	// the map which is created
	hashmap<int,hashset<int> > map;
	double radius;
	enum places_map_type { place_cell, place_place, cell_cell, cell_cell_centers };
	places_map_type pmt; // type of data saved in places_map
	
	// the data needed for it
		struct edge { // edge (between cells, used to store the Voronoi-polygons)
			int v1;
			int v2;
			int c1;
			int c2;
		};

		//0. cell coordinates and place coordinates
		hashmap<int,std::pair<double,double> > cell_coords;
		hashmap<int,std::pair<double,double> > place_coords;
		
		//1. list of vertex coordinates
		std::vector<double> vlon;
		std::vector<double> vlat;
		
		//2. vertex to cell matching
		hashmap<int,std::vector<int> > vertex_cell; //vertex to cell mapping (vertex i is included in the Voronoi polygons of cells in vertex_cell[i])
		
		//3. lines making up the Voronoi polygons (use these to test if the circles intersect with the lines)
		std::vector<edge> edges_cell;
		
		//4. cell to overlapping cell matching (only duplicates here)
		hashmap<int,std::vector<int> > cell_dup;
		hashmap<int,int> cell_dup_rev; //other way around
		
		// neighbors of cells (only among unique cell ids)
		hashmap<int,std::vector<int> > cell_neighb;
		
		
	// functions for reading the data
	int load_voronoi(FILE* f_vertex_coords, FILE* f_vertex_to_cell, FILE* f_edges, FILE* f_cell_unique);
	unsigned int load_cell_coords(FILE* cf, unsigned int header_skip, bool csvfile);
	unsigned int load_place_coords(FILE* cf, unsigned int header_skip, bool csvfile);
	
	// functions for creating the actual mapping
	unsigned int create_places_map(double radius);
	unsigned int create_places_map_self_cells(double radius, bool only_center);
	unsigned int create_places_map_self_places(double radius);
	
	// helper function for the previous ones (normally, this is only called internally)
	unsigned int create_places_map(double radius, const hashmap<int,std::pair<double,double> >& places,
		hashmap<int,hashset<int> >& places_map, bool add_closest);
	
	
	// print out the mapping
	void print_places_map(FILE* out) const {
		// file format: first line: type, second line: radius (both of these as "comments")
		// place ID, cell ID
		switch(pmt) {
			case place_cell:
				fprintf(out,"# place cell\n");
				break;
			case place_place:
				fprintf(out,"# place place\n");
				break;
			case cell_cell:
				fprintf(out,"# cell cell\n");
				break;
			case cell_cell_centers:
				fprintf(out,"# cell cell (centers)\n");
				break;
		}
		fprintf(out,"# %g\n",radius);
		const hashmap<int,hashset<int> >& places_map = map;
		for(hashmap<int,hashset<int> >::const_iterator it = places_map.begin(); it != places_map.end(); it++) {
			int pid = it->first;
			const hashset<int>& cells = it->second;
			for(hashset<int>::const_iterator it2 = cells.begin(); it2 != cells.end(); it2++)
				fprintf(out,"%d\t%d\n",pid,*it2);
		}
	}

	// helper functions for debugging
	int places_map_count(int place_id);
	int places_map_count(int place_id, int cell_id);
};


int main(int argc, char **argv)
{
	char* cellsfile = 0; //cell coordinates
	char* placesfile = 0; //place coordinates
	double radius = 500.0;
	
	//Voronoi tessellation
	//FILE* f_vertex_coords, FILE* f_vertex_to_cell, FILE* f_edges, FILE* f_cell_unique
	char* v_vertex_coords = 0;
	char* v_vertex_to_cell = 0;
	char* v_edges = 0;
	char* v_cell_unique = 0;
	
	bool cellsonly = false; //create the mapping only among the cells
	bool antennaonly = false; //as the previous, but only look at the antenna positions
	bool placesonly = false; //create the mapping only among the places
	
	//!! TODO: more rigorous processing of arguments, protect against invalid arguments
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'c':
			cellsfile = argv[i+1];
			break;
		case 'p':
			placesfile = argv[i+1];
			break;
		case 'd':
			radius = atof(argv[i+1]);
			break;
		case 'C':
			cellsonly = true;
			break;
		case 'A':
			antennaonly = true;
			cellsonly = true;
			break;
		case 'P':
			placesonly = true;
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
	
	if(placesonly) {
		if(placesfile == 0) {
			fprintf(stderr,"Error: missing input files!\n");
			return 1;
		}
	}
	else {
		if( (placesfile == 0 && cellsonly == false) || cellsfile == 0) {
			fprintf(stderr,"Error: missing input files!\n");
			return 1;
		}
		if( ! (v_vertex_coords && v_vertex_to_cell && v_edges && v_cell_unique) ) {
			fprintf(stderr,"Error: missing input files!\n");
			return 1;
		}
	}
	
	mapper cr;
	
	//read cell and place / stop coordinates, create matching between them
	FILE* f_vertex_coords = 0;
	FILE* f_vertex_to_cell = 0;
	FILE* f_edges = 0;
	FILE* f_cell_unique = 0;
	FILE* pf = 0;
	FILE* f_cell_coords = 0;
	
	if(!placesonly) {
		f_vertex_coords = fopen(v_vertex_coords,"r");
		f_vertex_to_cell = fopen(v_vertex_to_cell,"r");
		f_edges = fopen(v_edges,"r");
		f_cell_unique = fopen(v_cell_unique,"r");
		f_cell_coords = fopen(cellsfile,"r");
		if( ! ( f_vertex_coords && f_vertex_to_cell && f_edges && f_cell_unique && f_cell_coords ) ) {
			fprintf(stderr,"Error opening input files!\n");
			return 1;
		}
	}
	if(!cellsonly) {
		pf = fopen(placesfile,"r");
		if(!pf) {
			fprintf(stderr,"Error opening input files!\n");
			return 1;
		}
		if(cr.load_place_coords(pf,0,false) == 0) {
			fprintf(stderr,"Error reading place coordinates!\n");
			return 2;
		}
	}
	
	if(placesonly) {
		if(cr.create_places_map_self_places(radius) == 0) {
			fprintf(stderr,"Error creating place -- place matching!\n");
			return 2;
		}
	}
	else {
		if(cr.load_cell_coords(f_cell_coords,1,true) == 0) {
			fprintf(stderr,"Error reading cell coordinates!\n");
			return 2;
		}
		if(cr.load_voronoi(f_vertex_coords,f_vertex_to_cell,f_edges,f_cell_unique)) {
			fprintf(stderr,"Error reading Voronoi tessellation!\n");
			return 2;
		}
		
		if(cellsonly) {
			if(antennaonly) {
				if(cr.create_places_map_self_cells(radius,true) == 0) {
					fprintf(stderr,"Error creating antenna -- cell matching!\n");
					return 2;
				}
			}
			else if(cr.create_places_map_self_cells(radius,false) == 0) {
				fprintf(stderr,"Error creating cell -- cell matching!\n");
				return 2;
			}
		}
		else if(cr.create_places_map(radius) == 0) {
			fprintf(stderr,"Error creating place -- cell matching!\n");
			return 2;
		}
	}
	
	if(!placesonly) {
		fclose(f_vertex_coords);
		fclose(f_cell_coords);
		fclose(f_vertex_to_cell);
		fclose(f_edges);
		fclose(f_cell_unique);
	}
	if(!cellsonly) fclose(pf);
	
	cr.print_places_map(stdout);
	
	return 0;
}






// create place / stop -- cell mapping

// read cell coordinates from file, format should be id\tlon\tlat, skip the given amount of lines
// return value: number of lines read
static unsigned int read_coords(FILE* cf, hashmap<int,std::pair<double,double> >& coords, unsigned int header_skip, bool csvfile) {
	unsigned int r = 0;
	unsigned int line = 0;
	coords.clear();
	for(unsigned int j=0;j<header_skip;j++) {
		line++;
		int a;
		do {
			a = getc(cf);
		} while( ! (a == '\n' || a == EOF) );
		if(a == EOF) break;
	}
	if(feof(cf) || ferror(cf)) return r;
	
	while(1) {
		line++;
		int a = 0;
		int id = 0;
		double lon = 0.0;
		double lat = 0.0;
		if(csvfile) a = fscanf(cf,"%d,%lg,%lg",&id,&lon,&lat);
		else a = fscanf(cf,"%d%*[ \t]%lg%*[ \t]%lg",&id,&lon,&lat);
		if(a == EOF) break;
		if(a != 3) {
			fprintf(stderr,"crecords::load_cell_coords(): Invalid data in cell coordinate file, line %u!\n",line);
			break;
		}
		//!! note: we do not check if it was alrady there (i.e. multiple coordinates for the same cell in the same file)
		coords[id] = std::make_pair(lon,lat);
		r++;
		
		do {
			a = getc(cf);
		} while( ! (a == '\n' || a == EOF) );
		if(a == EOF) break;
	}
	
	return r;
}
		
unsigned int mapper::load_cell_coords(FILE* cf, unsigned int header_skip, bool csvfile) {
	return read_coords(cf, cell_coords, header_skip, csvfile);
}

unsigned int mapper::load_place_coords(FILE* cf, unsigned int header_skip, bool csvfile) {
	return read_coords(cf, place_coords, header_skip, csvfile);
}


// load pre-computed Voronoi tesselation (example: /data01/datasets/customer_matching)
// return value: 0: OK, >0: error (cell ID does not exists, or invalid input)
int mapper::load_voronoi(FILE* f_vertex_coords, FILE* f_vertex_to_cell, FILE* f_edges, FILE* f_cell_unique) {
	unsigned int line = 0;
	
	//clear all old data
	vlon.clear();
	vlat.clear();
	vertex_cell.clear();
	edges_cell.clear();
	cell_dup.clear();
	cell_dup_rev.clear();
	
	//1. load vertex coords
	//format: no header, tab or space separated
	while(1) {
		double lon = 0.0;
		double lat = 0.0;
		int a = fscanf(f_vertex_coords,"%lg%*[ \t]%lg",&lon,&lat);
		if(a == EOF) break;
		line++;
		if(a != 2) {
			fprintf(stderr,"Invalid data in vertex coord file line %u!\n",line);
			return 1;
		}
		vlon.push_back(lon);
		vlat.push_back(lat);
		
		do {
			a = getc(f_vertex_coords);
		} while( ! (a == '\n' || a == EOF) );
		if(a == EOF) break;
	}
	
	fprintf(stderr,"%lu vertices read\n",vlon.size());
	
	
	//2. load vertex to cell matching
	line = 0;
	while(1) {
		//first column: vertex_id, the rest are the cell ids (at least three)
		int vid = 0;
		int a = fscanf(f_vertex_to_cell,"%d",&vid);
		if(a == EOF) break;
		line++;
		if(a != 1) {
			fprintf(stderr,"Invalid data in the cell -- vertex matching file, line %d!\n",line);
			return 2;
		}
		
		if( vid >= vlon.size() || vid < -1 ) {
			fprintf(stderr,"Error: vertex id %d out of range (cell -- vertex matching file, line %d)!",vid,line);
			return 3;
		}
		
		if(vid > -1) { //ignore -1 (should not actually be in the file)
			std::vector<int>& vvec = vertex_cell[vid]; //note: new vector is created here
			if(vvec.size() > 0) { //already present, this is an error
				fprintf(stderr,"Error: vertex id %d appears multiple times in the cell -- vertex matching file (line %d)!\n",vid,line);
				return 4;
			}
			
			//read all cell ids
			while(1) {
				do {
					a = getc(f_vertex_to_cell);
				} while( a == ' ' || a == '\t' );
				if( a == '\n' || a == EOF ) break;
				ungetc(a,f_vertex_to_cell);
				int cid = 0;
				a = fscanf(f_vertex_to_cell,"%d",&cid);
				if(a != 1) {
					fprintf(stderr,"Invalid data in the cell -- vertex matching file, line %d!\n",line);
					return 2;
				}
				if(!cell_coords.count(cid)) {
					fprintf(stderr,"Invalid data in the cell -- vertex matching file, line %d: cell id %d not found!\n",line,cid);
					return 3;
				}
				vvec.push_back(cid); 
			}
			
			if(vvec.size() < 3) {
				fprintf(stderr,"Invalid data in cell -- vertex matching file, line %d: only %lu cells read for vertex id %d (expected >= 3)!\n",line,vvec.size(),vid);
				return 3;
			}
		} //vid > -1
		else {
			do {
				a = getc(f_vertex_to_cell);
			} while( ! (a == '\n' || a == EOF) );
		}
		
		//note: no need to advance to the end of the line, it is already reached in the previous loop
		if(a == EOF) break;
	}
	
	fprintf(stderr,"%u vertices read with matching cells\n",line);
	line = 0;
	
	//!! note: check if all vertices were read?
	
	//3. load line segments
	// f_edges file, each line contains four numbers: two cell ids, and two vertex ids
	while(1) {
		edge e;
		int a = fscanf(f_edges,"%d%*[ \t]%d%*[ \t]%d%*[ \t]%d",&e.c1,&e.c2,&e.v1,&e.v2);
		if(a == EOF) break;
		line++;
		if(a != 4) {
			fprintf(stderr,"Invalid data in the edges file, line %u!\n",line);
			return 6;
		}
		
		//check if cell and vertex ids are valid
		if( !cell_coords.count(e.c1) ) {
			fprintf(stderr,"Error in edges file, line %d: cell id %d not found!\n",line,e.c1);
			return 7;
		}
		if( !cell_coords.count(e.c2) ) {
			fprintf(stderr,"Error in edges file, line %d: cell id %d not found!\n",line,e.c2);
			return 7;
		}
		if(e.v1 >= (int)vlon.size() || e.v1 < -1) {
			fprintf(stderr,"Error in edges file, line %d: vertex id %d out of range!\n",line,e.v1);
			return 8;
		}
		if(e.v2 >= (int)vlon.size() || e.v2 < -1) {
			fprintf(stderr,"Error in edges file, line %d: vertex id %d out of range!\n",line,e.v2);
			return 8;
		}
		
		edges_cell.push_back(e);
		
		do {
			a = getc(f_edges);
		} while( ! (a == '\n' || a == EOF) );
		if(a == EOF) break;
	}
	
	fprintf(stderr,"%u edges read\n",line);
	line = 0;
	
	//4. load overlapping cells
	// f_cell_unique file, for each cell previous read (Voronoi polygons and edges) there is a line containing matching cells (including self match)
	// store only non-self matches
	while(1) {
		int id = 0;
		int id2 = 0;
		int a = fscanf(f_cell_unique,"%d%*[ \t]%d",&id,&id2);
		if(a == EOF) break;
		line++;
		if(a != 2) {
			fprintf(stderr,"Invalid data in overlapping cell matches file, line %d!\n",line);
			return 9;
		}
		if( !cell_coords.count(id) ) {
			fprintf(stderr,"Invalid data in overlapping cell matches file, line %d: cell id %d not found!\n",line,id);
			return 9;
		}
		
		/* if(id2 != id) */ cell_dup[id].push_back(id2);
		
		while(1) {
			do {
				a = getc(f_cell_unique);
			} while( a == ' ' || a == '\t' );
			if( a == '\n' || a == EOF ) break;
			ungetc(a,f_cell_unique);
			int cid = 0;
			a = fscanf(f_cell_unique,"%d",&cid);
			if(a != 1) {
				fprintf(stderr,"Invalid data in overlapping cell matches file, line %d!\n",line);
				return 9;
			}
			
			//!! note: check if coordinates are really the same?
			
			/* if(cid != id) */ cell_dup[id].push_back(cid); //note: this creates an empty vector if this id was not in there already
			// note: duplicates are not really a problem here, we create a hash set later
		}
		
		if(a == EOF) break;
	}
	
	//5. create reverse map
	for(hashmap<int, std::vector<int> >::const_iterator it = cell_dup.begin(); it != cell_dup.end(); ++it) {
		for(std::vector<int>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
			cell_dup_rev[*it2] = it->first;
		}
		cell_dup_rev[it->first] = it->first; //I'm not sure if this is necessary
	}
	
	
	//!! TODO: check if all cells have been read
	return 0;
}

// auxiliary functions for testing intersections

// calculate distance
// use spherical geometry

//!! TODO: CIRC = 2*M_PI*RADIUS ??
static const double RADIUS = 6371000.0;
static const double CIRC = 40000000.0;

static inline double dist_m(double lon1, double lat1, double lon2, double lat2) {
	lon1 = M_PI * lon1 / 180.0;
	lat1 = M_PI * lat1 / 180.0;
	lon2 = M_PI * lon2 / 180.0;
	lat2 = M_PI * lat2 / 180.0;

	double s1 = sin ((lat1 - lat2) / 2.0);
	double s2 = sin ((lon1 - lon2) / 2.0);
	double r1 = s1 * s1 + s2 * s2 * cos (lat1) * cos (lat2);
	return RADIUS * 2.0 * asin (sqrt (r1));
}

// check if a circle and a line segment intersect
// note: use a local Euclidean geometry, scale longitudes by 1/cos(lat) factor
// will be precise only on small scales
// clon, clat: circle center coordinates
// r: circle radius (meters)
// elon1, elat1, elon2, elat2: ends of the line segment
static inline bool circle_edge_intersect(double clon, double clat, double r, double elon1, double elat1, double elon2, double elat2) {
	if(r >= CIRC/2.0) return false;
	double sf = cos( clat * M_PI / 180.0 );
	if(sf < 0.001) {
		fprintf(stderr,"circle_edge_intersect(): error: circle is too close to poles!\n");
		return false;
	}
	
	clon /= sf;
	elon1 /= sf;
	elon2 /= sf;
	
	r = 360.0 * r / CIRC;
	
	//note: coordinates of the line segment:
	// x = elon1 + t*(elon2 - elon1)
	// y = elat1 + t*(elat2 - elat1)
	//where t \in [0,1]
	// equation defining intersections:
	// ( x - clon )^2 + ( y - clat )^2 = r^2
	//we have an intersection if this can be solved and the solution t is between 0 and 1
	
	double a = (elon1 - elon2)*(elon1 - elon2) + (elat1 - elat2)*(elat1 - elat2);
	double b = 2*(elon2 - elon1)*(elon1 - clon) + 2*(elat2 - elat1)*(elat1 - clat);
	double c = (elon1 - clon)*(elon1 - clon) + (elat1 - clat)*(elat1 - clat) - r*r;
	
	double D = b*b - 4*a*c;
	if(D < 0) return false;
	
	double t1 = ( - b + sqrt(D) ) / (2*a);
	if(t1 >= 0.0 && t1 <= 1.0) return true;
	double t2 = ( - b - sqrt(D) ) / (2*a);
	if(t2 >= 0.0 && t2 <= 1.0) return true;
	return false;
}

// check if a circle and a half-line intersect
// note: use a local Euclidean geometry, scale longitudes by 1/cos(lat) factor
// will be precise only on small scales
// clon, clat: circle center coordinates
// r: circle radius (meters)
// elon1, elat1: start of the line
// elon2, elat2: any point on the line (line goes in this direction from start); practically, the point halfway between the two cell coordinates
static inline bool circle_halfline_intersect(double clon, double clat, double r, double elon1, double elat1, double elon2, double elat2) {
	if(r >= CIRC / 2.0) return false;
	double sf = cos( clat * M_PI / 180.0 );
	if(sf < 0.001) {
		fprintf(stderr,"circle_edge_intersect(): error: circle is too close to poles!\n");
		return false;
	}
	
	clon /= sf;
	elon1 /= sf;
	elon2 /= sf;
	
	r = 360.0 * r / CIRC;
	
	//note: coordinates of the half-line:
	// x = elon1 + t*(elon2 - elon1)
	// y = elat1 + t*(elat2 - elat1)
	//where t >= 0
	// equation defining intersections:
	// ( x - clon )^2 + ( y - clat )^2 = r^2
	//we have an intersection if this can be solved and the solution is t >= 0
	double a = (elon1 - elon2)*(elon1 - elon2) + (elat1 - elat2)*(elat1 - elat2);
	double b = 2*(elon2 - elon1)*(elon1 - clon) + 2*(elat2 - elat1)*(elat1 - clat);
	double c = (elon1 - clon)*(elon1 - clon) + (elat1 - clat)*(elat1 - clat) - r*r;
	
	double D = b*b - 4*a*c;
	if(D < 0) return false;
	
	double t1 = ( - b + sqrt(D) ) / (2*a);
	if(t1 >= 0.0) return true;
	double t2 = ( - b - sqrt(D) ) / (2*a);
	if(t2 >= 0.0) return true;
	return false;
}


// get direction (coordinates of a point on the line) of a halfline starting at v between cells p1 and p2
// c is the coordinates of the center, sf is the scaling factor to use for the longitudes (sf = cos(c.second * M_PI / 180.0))
static inline std::pair<double,double> halfline_get_coords(const std::pair<double,double>& v, const std::pair<double,double>& p1,
		const std::pair<double,double>& p2,	const std::pair<double,double>& c, double sf) {
	double elon2 = (p1.first + p2.first) / 2.0;
	double elat2 = (p1.second + p2.second) / 2.0;
				
	//test if (elon2, elat2) is closer to the center than (vlon[v1],vlat[v1])
	double dx1 = elon2 - v.first;
	double dy1 = elat2 - v.second;
	double dx2 = c.first - v.first;
	double dy2 = c.second - v.second;
	
	// scale longitudes in a local Euclidean geometry -- this only affects the calculation of the scalar product; all other calculations here
	//	are linear (i.e. additions, subtractions and scalings), which are not affected by an arbitrary scaling factor
	if(dx1*dx2 / (sf*sf) + dy1*dy2 > 0) { //these point in the same direction, edge coordinates need to be swapped
		elon2 = v.first - dx1;
		elat2 = v.second - dy1;
	}
	
	return std::make_pair(elon2,elat2);
}


// check if two lines intersect
// use Euclidean geometry WITHOUT rescaling coordinates (it should already be scaled previously)
// the lines are defined by segments e1 -> e2 and f1 -> f2
// returns true if the intersection is between e1 and e2 and f1 and f2 (on both lines), or between f1 and f2 and in the direction of
//		e2 from e1 if ehalfline == true
// false otherwise (outside or parallel lines)
static inline int lines_intersect(double e1lon, double e1lat, double e2lon, double e2lat, double f1lon, double f1lat,
		double f2lon, double f2lat, int ehalfline, int fhalfline) {
	
	// equation defining the instersection: e1 + (e2-e1)t = f1 + (f2-f1)s
	// t and s are the scalar unknowns, each eq. has two components
	// the equation has a solution if the determinant is nonzero
	double vlon = e2lon - e1lon;
	double vlat = e2lat - e1lat;
	double ulon = f2lon - f1lon;
	double ulat = f2lat - f1lat;
	double alon = f1lon - e1lon;
	double alat = f1lat - e1lat;
	// equation with these variables:
	// v*t - u*s = a
	// solution:
	// D = u_1 * v_2 - v_1 * u_2
	// t = ( u_1 * a_2 - u_2 * a_1 ) / D
	// s = ( v_1 * a_2 - v_2 * a_1 ) / D
	double d = ulon * vlat - vlon * ulat;
	if(d < 1e-20) return 0; // (nearly) parallel lines
	double t = ( ulon * alat - alon * ulat ) / d;
	double s = ( vlon * alat - alon * vlat ) / d;

	if(s < 0 || t < 0) return 0;
	if(s > 1 || !fhalfline) return 0;
	if(t > 1 || !ehalfline) return 0;
	return 1;
}

// check if point p is closer than r to the line segment e1->e2 (or halfline from e1 towards e2 if halfline == true)
//	the point needs to fall between e1 and e2 on the line
// use Euclidean geometry WITHOUT rescaling coordinates (it should already be scaled previously)
static inline int point_line_check(double e1lon, double e1lat, double e2lon, double e2lat, double plon, double plat, double r, int halfline) {
	// p1 = p-e1
	// v = e2-e1
	double p1lon = plon - e1lon;
	double p1lat = plat - e1lat;
	double vlon = e2lon - e1lon;
	double vlat = e2lat - e1lat;
	// distance from the line:
	// d^2 = p1^2 - (p1*v)^2 / v^2
	// "shadow" on the line:
	// t = (p1*v) / |v|
	double p1v = p1lon * vlon + p1lat * vlat;
	if(p1v < 0) return 0; // wrong direction
	double v2 = vlon * vlon + vlat * vlat;
	if(halfline == false) if(p1v > v2) return 0; // outside
	double p12 = p1lon * p1lon + p1lat * p1lat;
	double d2 = p12 - p1v*p1v / v2;
	if(d2 > r*r) return 0;
	return 1;
}


static int places_map_debug_c1 = -1;
static int places_map_debug_c2 = -1;

static unsigned int places_map_insert(hashmap<int,hashset<int> >& places_map, int c1, int c2, bool symm) {
	unsigned int r = 0;
	if( (c1 == places_map_debug_c1 && c2 == places_map_debug_c2) || (c2 == places_map_debug_c1 && c1 == places_map_debug_c2 && symm) ) {
		fprintf(stderr,"debug\n"); // set a breakpoint here and set the debug variables
	}
	if(places_map[c1].insert(c2).second) r++;
	if(symm) if(places_map[c2].insert(c1).second) r++;
	return r;
}
		



unsigned int mapper::create_places_map(double radius, const hashmap<int,std::pair<double,double> >& places,
		hashmap<int,hashset<int> >& places_map, bool add_closest) {
	//~ unsigned int n = places.size();
	if(radius <= 0.0) return 0;
	
	unsigned int r = 0;
	
	places_map.clear(); // make sure it's empty
	
	//calculate the center (mean of all cell coordinates) to use for half-lines
	double clon = 0.0;
	double clat = 0.0;
	unsigned int nc = 0;
	for(hashmap<int,std::pair<double,double> >::const_iterator it = cell_coords.begin(); it != cell_coords.end(); ++it) {
		clon += it->second.first;
		clat += it->second.second;
		nc++;
	}
	clon /= (double)nc;
	clat /= (double)nc;
	double sf = cos(clat * M_PI / 180.0);
	
	for(hashmap<int,std::pair<double,double> >::const_iterator it = places.begin(); it != places.end(); ++it) {
		int pid = it->first;
		double lon = it->second.first;
		double lat = it->second.second;
		
		if(add_closest) {
			// 0. check all cells, add the closest one (the Voronoi cell containing the place)
			double dmin = 40000000.0;
			int cmin = -1;
			for(hashmap<int,std::pair<double,double> >::const_iterator it = cell_coords.begin(); it != cell_coords.end(); ++it) {
				double d1 = dist_m(lon, lat, it->second.first, it->second.second);
				if(d1 < dmin) {
					dmin = d1;
					cmin = it->first;
				}
			}
			if(cmin == -1) {
				//error, this should never happen, distances on earth are < 20000km
				fprintf(stderr,"crecords::create_places_map(): error in calculating distances!\n");
				return 0;
			}
			//create map for this place, add the closest
			hashmap<int,int>::const_iterator it2 = cell_dup_rev.find(cmin);
			if(it2 != cell_dup_rev.end()) {
				cmin = it2->second;
				hashmap<int, std::vector<int> >::const_iterator it3 = cell_dup.find(cmin);
				if(it3 != cell_dup.end()) {
					const std::vector<int>& c2 = it3->second;
					for(std::vector<int>::const_iterator itc2 = c2.begin(); itc2 != c2.end(); ++itc2)
							r += places_map_insert(places_map,pid,*itc2,false);
				}
			}
			r += places_map_insert(places_map,pid,cmin,false);
		}
		
		
		// 1. check all vertices
		for(unsigned int j=0;j<vlon.size();j++) {
			double d1 = dist_m(lon,lat,vlon[j],vlat[j]);
			if(d1 < radius) {
				const std::vector<int>& c1 = vertex_cell[j];
				for(unsigned int k=0;k<c1.size();k++) {
					int cid = c1[k];
					//insert cell id
					r += places_map_insert(places_map,pid,cid,false); //note: if it is already there, insert() is a no-op
					//check for duplicate cell_ids and insert them too
					hashmap<int, std::vector<int> >::const_iterator it = cell_dup.find(cid);
					if(it != cell_dup.end()) {
						const std::vector<int>& c2 = it->second;
						for(std::vector<int>::const_iterator itc2 = c2.begin(); itc2 != c2.end(); ++itc2)
							r += places_map_insert(places_map,pid,*itc2,false);
					}
				}
			} // d1 < radius (insert cells)
		}
		
		// 2. check all edges
		for(unsigned int j=0;j<edges_cell.size();j++) {
			int c1 = edges_cell[j].c1;
			int c2 = edges_cell[j].c2;
			int v1 = edges_cell[j].v1;
			int v2 = edges_cell[j].v2;
			bool intersects;
			if(v1 == -1 || v2 == -1) {
				// this is a half-line
				if(v1 == -1) { // make v2 the unknown
					int tmp = v1;
					v1 = v2;
					v2 = tmp;
				}
				// calculate the direction of the line
				const std::pair<double,double> e2 = halfline_get_coords(std::make_pair(vlon[v1],vlat[v1]),cell_coords[c1],cell_coords[c2],
					std::make_pair(clon,clat),sf);
				intersects = circle_halfline_intersect(lon, lat, radius, vlon[v1], vlat[v1], e2.first, e2.second);
			}
			else intersects = circle_edge_intersect(lon, lat, radius, vlon[v1], vlat[v1], vlon[v2], vlat[v2]);
			if(intersects) {
				// add both cells
				r += places_map_insert(places_map,pid,c1,false);
				r += places_map_insert(places_map,pid,c2,false);
				//check for duplicate cell_ids and insert them too
				hashmap<int, std::vector<int> >::const_iterator it = cell_dup.find(c1);
				if(it != cell_dup.end()) {
					const std::vector<int>& cv1 = it->second;
					for(std::vector<int>::const_iterator itcv1 = cv1.begin(); itcv1 != cv1.end(); ++itcv1)
						r += places_map_insert(places_map,pid,*itcv1,false);
				}
				hashmap<int, std::vector<int> >::const_iterator it2 = cell_dup.find(c2);
				if(it2 != cell_dup.end()) {
					const std::vector<int>& cv2 = it2->second;
					for(std::vector<int>::const_iterator itcv2 = cv2.begin(); itcv2 != cv2.end(); ++itcv2)
						r += places_map_insert(places_map,pid,*itcv2,false);
				}
			}
		}
	}
	
	this->radius = radius;
	return r;
}


// auxiliary functions for debugging
int mapper::places_map_count(int place_id) {
	return map.count(place_id);
}
int mapper::places_map_count(int place_id, int cell_id) {
	hashmap<int,hashset<int> >::const_iterator it = map.find(place_id);
	if(it == map.end()) return 0;
	else return it->second.count(cell_id);
}


// if searchreverse == true, the mappings are created reversed (cell_id -> place_id), this should be used if this class contains transportation records
unsigned int mapper::create_places_map(double radius) {
	pmt = place_cell;
	this->radius = radius;
	return create_places_map(radius, place_coords, map, true);
}
	


unsigned int mapper::create_places_map_self_places(double radius) {
	unsigned int r = 0;
	
	// calculate the distance between all pairs
	for(hashmap<int,std::pair<double,double> >::const_iterator it = place_coords.begin(); it != place_coords.end(); ++it) {
		int pid = it->first;
		
		// add self first
		r += places_map_insert(map,pid,pid,false);
		
		bool found = false;
		for(hashmap<int,std::pair<double,double> >::const_iterator it2 = place_coords.begin(); it2 != place_coords.end(); ++it2) {
			int pid2 = it2->first;
			if(pid2 >= pid) continue;
			
			double r1 = dist_m( it->second.first, it->second.second, it2->second.first, it2->second.second );
			if(r1 <= radius) r += places_map_insert(map,pid,pid2,true);
		}
	}
	
	pmt = place_place;
	this->radius = radius;
	return r;
}



unsigned int mapper::create_places_map_self_cells(double radius, bool only_center) {

	// two modes of operation:
	//	1. only_center == true
	//		in this case, the cell coordinates (from cell_coords) are passed to create_places_map() with the original hashmap
	//		(i.e. the cells in the given radius of the antennas are searched)
	//		duplicate cells are handled after running and added as duplicates as well
	//	2. only_center == false
	//		in this case, the vertex coordinates are passed to create_places_map() which creates along with a temporary hashmap which is
	//		filled with which vertex is close to which cell
	//		the original hashmap is then filled based on the vertex to cell matching
	
	unsigned int r = 0;
	hashmap<int,hashset<int> >& places_map = map; // alias for the real map created
	if(only_center) {
		// 1. search for neighbors of the cells
		r += create_places_map(radius,cell_coords,map,false); // last parameter: do not add the closest cell, only cells which are indeed in the radius
		
		// 2. make the places_map symmetric
		std::vector<std::pair<int,int> > to_insert; // we store missing matches here (I'm not sure if inserting to a hashmap while iterating is a good idea)
		for(hashmap<int,hashset<int> >::const_iterator it = map.begin(); it != map.end(); ++it) {
			const hashset<int>& values = it->second;
			int cid1 = it->first;
			for(hashset<int>::const_iterator it2 = values.begin(); it2 != values.end(); ++it2) {
				int cid2 = *it2;
				hashmap<int,hashset<int> >::const_iterator it3 = map.find(cid2);
				if(it3 == map.end()) to_insert.push_back(std::make_pair(cid2,cid1)); //surely this is missing
				else if(it3->second.count(cid1) == 0) to_insert.push_back(std::make_pair(cid2,cid1)); //check if the set of cid2 contains cid1
			}
		}
		
		for(std::vector<std::pair<int,int> >::const_iterator it = to_insert.begin(); it != to_insert.end(); ++it)
			r += places_map_insert(map,it->first,it->second,false);
		
		pmt = cell_cell_centers;
	}
	else {
		// 1. create a hashmap numbering the vertices (needed for create_places_map()) and an empty temporary places_map
		hashmap<int,std::pair<double,double> > vertex_coords_map;
		for(unsigned int i=0;i<vlon.size();i++) vertex_coords_map[(int)i] = std::make_pair(vlon[i],vlat[i]);
		hashmap<int,hashset<int> > places_map_tmp;
		
		// 2. search for vertex neighbors
		create_places_map(radius,vertex_coords_map,places_map_tmp,false); // last parameter: do not add the closest cell, only cells which are indeed in the radius
		
		
		// 3. create the real places_map based on vertex to cell mapping
		for(hashmap<int,hashset<int> >::const_iterator it = places_map_tmp.begin(); it != places_map_tmp.end(); ++it) {
			int vid = it->first;
			hashmap<int,std::vector<int> >::const_iterator it2 = vertex_cell.find(vid);
			if(it2 != vertex_cell.end())
				for(std::vector<int>::const_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
					for(hashset<int>::const_iterator it4 = it->second.begin(); it4 != it->second.end(); ++it4)
						// *it3 and *it4 are the two cell IDs which are a match here
						r += places_map_insert(places_map,*it3,*it4,true); // note: keep track of the matches found
					
		}
		
		// 4. search for the neighbors of all line segments
		// calculate the center (mean of all cell coordinates) to use for half-lines and scaling factors
		// in this part, we only work with scaled coordinates
		double clon = 0.0;
		double clat = 0.0;
		unsigned int nc = 0;
		for(hashmap<int,std::pair<double,double> >::const_iterator it = cell_coords.begin(); it != cell_coords.end(); ++it) {
			clon += it->second.first;
			clat += it->second.second;
			nc++;
		}
		clon /= (double)nc;
		clat /= (double)nc;
		double sf = cos(clat * M_PI / 180.0);
		double r1 = 360.0 * radius / CIRC; // rescale distance to scaled degrees
		for(uint i=0;i<edges_cell.size();i++) {
			int v1 = edges_cell[i].v1;
			int v2 = edges_cell[i].v2;
			int c1 = edges_cell[i].c1;
			int c2 = edges_cell[i].c2;
			
			// 4.1 add c1 and c2 to each other's list (neighbors always match in this case)
			//	(possibly this was already done in the previous step)
			r += places_map_insert(places_map,c1,c2,true);
			
			double v1lon,v1lat,v2lon,v2lat;
			int halfline = 0;
			if(v1 == -1 || v2 == -1) {
				halfline = 1;
				if(v1 == -1) {
					v1 = v2;
					v2 = -1;
				}
				v1lon = vlon[v1] / sf;
				v1lat = vlat[v1];
				//~ std::pair<double,double> v2p = halfline_get_coords(std::make_pair(v1lon,v1lat),cell_coords.at(c1),cell_coords.at(c2),
					//~ std::make_pair(clon,clat),sf); //!!
				std::pair<double,double> v2p = halfline_get_coords(std::make_pair(vlon[v1],v1lat),cell_coords.at(c1),cell_coords.at(c2),
					std::make_pair(clon,clat),sf);
				v2lon = v2p.first / sf;
				v2lat = v2p.second;
			}
			else {
				v1lon = vlon[v1] / sf;
				v1lat = vlat[v1];
				v2lon = vlon[v2] / sf;
				v2lat = vlat[v2];
			}
			
			// check all points to be close
			for(size_t j=0;j<vlon.size();j++) if( ! (j == v1 || j == v2 ) ) {
				double v3lon = vlon[j] / sf;
				double v3lat = vlat[j];
				if(point_line_check(v1lon,v1lat,v2lon,v2lat,v3lon,v3lat,r1,halfline)) {
					// match, add all candidates
					std::vector<int>& cells1 = vertex_cell.at((int)j); // exception if not found
					for(std::vector<int>::const_iterator it = cells1.begin(); it != cells1.end(); ++it) {
						r += places_map_insert(places_map,c1,*it,true);
						r += places_map_insert(places_map,c2,*it,true);
					}
				}
			}
			
			// create parallel line segments, check for intersection with all other edges
			double dx = v2lat - v1lat; // vector at right angles to the line segment
			double dy = v1lon - v2lon;
			// normalize length
			double dr = r1 / sqrt(dx*dx + dy*dy);
			dx *= dr;
			dy *= dr;
			
			double f1lon = v1lon + dx;
			double f1lat = v1lat + dy;
			double f2lon = v2lon + dx;
			double f2lat = v2lat + dy;
			double g1lon = v1lon - dx;
			double g1lat = v1lat - dy;
			double g2lon = v2lon - dx;
			double g2lat = v2lat - dy;
			
			// check all other lines for intersection
			for(uint k=i+1;k<edges_cell.size();k++) {
				int halfline2 = 0;
				if(edges_cell[k].v1 == -1 || edges_cell[k].v2 == -1) halfline2 = 1;
				double u1lon,u2lon,u1lat,u2lat;
				if(edges_cell[k].v1 == -1) {
					u1lon = vlon[ edges_cell[k].v2 ] / sf;
					u1lat = vlat[ edges_cell[k].v2 ];
				}
				else {
					u1lon = vlon[ edges_cell[k].v1 ] / sf;
					u1lat = vlat[ edges_cell[k].v1 ];
				}
				if(edges_cell[k].v2 == -1) {
					std::pair<double,double> tmp = halfline_get_coords(std::make_pair(u1lon,u1lat), cell_coords.at( edges_cell[k].c1 ),
						cell_coords.at( edges_cell[k].c2 ), std::make_pair(clon,clat), sf);
					u2lon = tmp.first / sf;
					u2lat = tmp.second;
				}
				else {
					u2lon = vlon[ edges_cell[k].v2 ] / sf;
					u2lat = vlat[ edges_cell[k].v2 ];
				}
				
				// check intersection
				bool match = false;
				if(lines_intersect(u1lon,u1lat,u2lon,u2lat,f1lon,f1lat,f2lon,f2lat,halfline,halfline2)) match = true;
				else if (lines_intersect(u1lon,u1lat,u2lon,u2lat,g1lon,g1lat,g2lon,g2lat,halfline,halfline2)) match = true;
				if(match) {
					int c3 = edges_cell[k].c1;
					int c4 = edges_cell[k].c2;
					// four matches: c1 -- c3, c1 -- c4, c2 -- c3, c2 -- c4 (also c1 -- c2 and c3 -- c4, but those are already taken care of)
					r += places_map_insert(places_map,c1,c3,true) + places_map_insert(places_map,c1,c4,true) +
							places_map_insert(places_map,c2,c3,true) + places_map_insert(places_map,c2,c4,true);
				}
			}
		}
		
		pmt = cell_cell;
	}
	
	// 4. / 5. make sure that every cell is matched with itself (only non-duplicate cells, duplicates will be added later on)
	for(hashmap<int,std::vector<int> >::const_iterator it = cell_dup.begin(); it != cell_dup.end(); ++it) {
		hashset<int>& values = places_map[it->first];
		if(values.insert(it->first).second) r++;
	}
	
	// 5. / 6. create copies for the duplicates
	// 5.1. for every match look up possible duplicates (including duplicats of the cell itself)
	std::vector<int> toadd; //vector to hold temporary ids to add to a hashset
	for(hashmap<int,hashset<int> >::iterator it = places_map.begin(); it != places_map.end(); ++it) {
		int c1 = it->first;
		hashset<int>& values = it->second;
		if(values.insert(c1).second) r++;
		for(hashset<int>::const_iterator it3 = values.begin(); it3 != values.end(); ++it3) {
			hashmap<int,std::vector<int> >::const_iterator it2 = cell_dup.find(*it3);
			if(it2 != cell_dup.end()) {
				toadd.insert(toadd.end(),it2->second.begin(),it2->second.end());
			}
		}
		for(std::vector<int>::const_iterator it4 = toadd.begin(); it4 != toadd.end(); ++it4) {
			if(values.insert(*it4).second) r++;
		}
		toadd.clear();
	}
	
	// 5.2 make copies of the sets for each duplicate
	for(hashmap<int,std::vector<int> >::const_iterator it = cell_dup.begin(); it != cell_dup.end(); ++it) {
		int cid = it->first;
		
		for(std::vector<int>::const_iterator it3 = it->second.begin(); it3 != it->second.end(); ++it3) {
			places_map[*it3] = places_map[cid];
			r += places_map[cid].size();
		}
	}
	
	this->radius = radius;
	return r;
}


