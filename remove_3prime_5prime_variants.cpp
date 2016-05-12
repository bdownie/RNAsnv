#include <stdio.h>  
#include "api/BamReader.h"

#include <getopt.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>

#include <iomanip> 

//#include <boost/algorithm/string.hpp>


using namespace BamTools; 
using namespace std;
      
int main(int argc, char *argv[])  
{  
	char c;
	extern char *optarg;

	//int clip_dist = -1;
	unsigned int min_alt_count = 1;
	//unsigned int min_splice_count = 1;

	char *VCF_FILE = NULL;
	char *BAM_FILE = NULL;
	//char *REF_FILE = NULL;

	int minimum_median_end_dist = -1;
	unsigned int min_homopolymer_length = 5;

	setprecision(3);

	while ((c = getopt (argc, argv, "m:b:v:r:d:h")) != -1) {
		switch(c) {
			case 'b':
				BAM_FILE = optarg;
				break;
			case 'v':
				VCF_FILE = optarg;
				break;
	//		case 'r':
	//			REF_FILE = optarg;
	//			break;
			case 'd':
				minimum_median_end_dist = atoi(optarg);
				break;
			case 'm':
				min_homopolymer_length = atoi(optarg);
				break;
			case 'h':
				cerr << "Usage: -b <bam> -v <vcf> -r <ref> -d <min median end distance>" << endl;
				return -1;
			default:
				cerr << "Usage: -b <bam> -v <vcf> -r <ref> -d <min median end distance>" << endl;
				return -1;
		}
	}

	if (NULL == BAM_FILE) { 
		cerr << "Usage: -b <bam> -v <vcf> -r <ref> -d <min median end distance>" << endl;
		return -1;
	}

	istream *in;
	ifstream v;
	if (NULL == VCF_FILE) { 
		in = &cin;
	}
	else {
		v.open(VCF_FILE);
		in = &v;
	}


	string line;
	stringstream tmp_ss;


	BamReader b;
	b.Open(BAM_FILE);
	//BamIndex::IndexType idx = BamIndex::BAMTOOLS;
	//if ((!b.LocateIndex(BamIndex::BAMTOOLS)) && (!b.CreateIndex(BamIndex::BAMTOOLS))) {
		if ((!b.LocateIndex(BamIndex::STANDARD)) && (!b.CreateIndex(BamIndex::STANDARD))) {
			cerr << "Error: " << b.GetErrorString() << endl;
			return 1;
		}
	//}
	//string nm_string("NM");
	while (getline(*in,line)) {
		if (line[0] == '#') { cout << line << endl; continue; }
		//vector <string> split_line;
		//boost::split(split_line,line,boost::is_any_of("\t "));
		vector <string> split_line;
		stringstream tmp_ss;
		tmp_ss.str(line);
                unsigned int fieldnum = 0;
		while (tmp_ss.good() && (fieldnum < 8)) { 
			string tmp;
			tmp_ss >> tmp;
			split_line.push_back(tmp);
		}
		string var = split_line[4];
                if (var == ".") {
                    if (split_line[7].find("END=") != string::npos) { 
                        cout << line << endl;
                    }
                    continue;
                }

		string chr = split_line[0];
		string site = split_line[1];
		string trash;


		//tmp_ss >> chr;
		//tmp_ss >> site;
		//tmp_ss >> trash;
		//tmp_ss >> trash;
		//tmp_ss >> var;

		int isite = atoi(site.c_str());
		isite--;

		int ichr = b.GetReferenceID(chr);
		BamRegion br(ichr,isite,ichr,isite);
		if (!b.SetRegion(br)) { 
			cerr << "Couldn't jump to " << ichr << " " << isite << endl;
			cerr << b.GetErrorString() << endl;
			return 2;
			//RefVector r = b.GetReferenceData();
		}

		BamAlignment myalign;
		unsigned int alt_count = 0;
		//unsigned int has_n = 0;
		//unsigned int only_indel = 0;

		vector <float> median_edit_dist;
		vector <unsigned int> MedianReadEndDist;
		//unsigned int MaxDistToReadEnd = 0;

		while(b.GetNextAlignment(myalign)) {
			int32_t edit_dist = -1;
			float edit_dist_normalized = 0;
			if	(myalign.HasTag("NM") && myalign.GetTag("NM",edit_dist)) {
				if (edit_dist == 0) { continue; }
				edit_dist_normalized = (float)edit_dist/myalign.Length;
			}
			else if (myalign.HasTag("nM") && myalign.GetTag("nM",edit_dist)) {
				if (edit_dist == 0) { continue; }
				edit_dist_normalized = (float)edit_dist/myalign.Length;
				if (myalign.IsMateMapped()) {
					edit_dist_normalized = edit_dist_normalized/2;
				}
			}

			if (minimum_median_end_dist < 0) { 
				minimum_median_end_dist = myalign.Length/10;
				minimum_median_end_dist--;
				cerr << "Setting min median end distance to " << myalign.Length << " " << minimum_median_end_dist << endl;
			}
			vector <int> soft_clips;
			vector <int> soft_read_loc;
			vector <int> soft_genome_loc;

			string bases = myalign.AlignedBases;

			unsigned int loc = isite - myalign.Position;


			if (bases[loc] == 'N') { 
				//has_n++;
				continue;
			}
			unsigned int loc2 = myalign.GetEndPosition() - isite;
			if (loc < 0 || loc2 < 0) { 
				cerr << "Bad position for read: " << myalign.Name << endl;
				return 2;
			}

			vector<CigarOp> cigar = myalign.CigarData;
			bool skip_read = 0;
			unsigned int n_len = 0;
			for (vector<CigarOp>::iterator i=cigar.begin(); i != cigar.end(); i++) {
				if ((i->Type == 'I') || (i->Type == 'D')) {
					//only_indel++;
					skip_read = 1;
					continue;
				}
				if (i->Type == 'N') {
					n_len += i->Length;
				}
			}

			// Iterator backward for homopolymer search if on forward strand
			int strand_iterator = -1;
			if (myalign.IsReverseStrand()) {
				strand_iterator = 1;
			}

			unsigned int num_test = 0;
			bool homopolymer_found = 1;
			if (min_homopolymer_length > 0) { 
				for (int homopolymer_i = loc + strand_iterator; num_test < min_homopolymer_length; homopolymer_i += strand_iterator) {

					if ((homopolymer_i < 0) || (homopolymer_i > myalign.Length) || bases[homopolymer_i] != var[0]) { 
						
						
						homopolymer_found = 0; num_test = min_homopolymer_length + 1; }
				}
				if (homopolymer_found) { 
					
					skip_read = 1; }
			}

			if (!skip_read && (bases[loc] == var[0])) {
				median_edit_dist.push_back(edit_dist_normalized);

				if (loc > n_len) { loc -= n_len; }
				if (loc2 > n_len) { loc2 -= n_len; }

				if (loc2 < loc) { loc = loc2; }
				MedianReadEndDist.push_back(loc);
			}


			alt_count++;



		}

		


		if (alt_count >= min_alt_count && (median_edit_dist.size() > min_alt_count)) {
			sort(median_edit_dist.begin(), median_edit_dist.end());
			sort(MedianReadEndDist.begin(), MedianReadEndDist.end());
			float median = 0;
			unsigned int size = median_edit_dist.size();
			//cout << size << endl;
			if (size % 2 == 0) {
    			median = (median_edit_dist[size / 2 - 1] + median_edit_dist[size / 2]) / 2;
  			}
			else { 
      			median = median_edit_dist[size / 2];
			}
			median_edit_dist.clear();
			split_line[7].append(";MedianEditDist=");
			stringstream tmp;
			tmp.precision(3);
			tmp << median;


			split_line[7].append(tmp.str());
//			cout << size << " " << median<< " " << split_line[7] << endl;

			median = 0;
			size = MedianReadEndDist.size();
			if (size % 2 == 0) {
    			median = (MedianReadEndDist[size / 2 - 1] + MedianReadEndDist[size / 2]) / 2;
  			}
			else { 
      			median = MedianReadEndDist[size / 2];
			}
			//line = boost::join(split_line,'\t');
			if (median >= minimum_median_end_dist) {  
				split_line[7].append(";MedianReadEndDist=");
				tmp.str("");
				tmp << median;
				split_line[7].append(tmp.str());
				split_line[7].append(";MaxReadEndDist=");
				tmp.str("");
				tmp << MedianReadEndDist[MedianReadEndDist.size() - 1];
				MedianReadEndDist.clear();
				split_line[7].append(tmp.str());

				for (unsigned int it = 0; it < split_line.size(); it++) { 
					if (it != 0) { cout << "\t"; }
					cout << split_line[it];
				}
				cout << endl;
			}
		}
	}
    return 0;  
}  

