/*
 *  fq2sc.cpp
 *  
 *
 *  Created by Jianhong Ou on 3/30/12.
 *  Copyright 2012 UMASSMED. All rights reserved.
 *
 */

/*
 * convert fastq file to sequence count file
 * s2c file
 * first column is sequence
 * second column is the count
 * sample:
 *		AATTGGCC	10
 *		GGCCTTAA	8
 */

#include <zlib.h> /* for gziped file */
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <unistd.h> /* for getopt */
#include "kseq.h" /* for parse fastq,fasta */
#include "fq2sc.h" /* for parse fastq to s2c */

using namespace std;

void usage(char* cmd){
	fprintf(stderr,"Usage: %s -i input -o output\n"
			"\ti\t<fastq file name>\n"
			"\to\t<output file name>\n"
			"\tq\t<if is set, do QA filter>\n"
			"\tc\t<quality cutoff threthold, need -q>\n", cmd);
}

// declare the type of file handler and the read() function 
KSEQ_INIT(gzFile, gzread)

// calculate qulity for each position
void getQA(string s, map<int,double>&qa, int l, int &max, int &min){
	map<int,double>::iterator iit;
	for (int i=0; i<s.length(); i++) {
		iit=qa.find(i);
		if (iit!=qa.end()) {
			qa[i]=(qa[i]*(l-1)+s[i])/l;
		}else {
			qa[i]=s[i];
		}
		//calculate the encodeing of qualities. Phred+33 or Phred/Solexa + 64?
		//so if max > 74 and min > 58 : +64; otherwise +33;
		if (min>(int)s[i]) {
			min=(int)s[i];
		}
		if (max<(int)s[i]) {
			max=(int)s[i];
		}
	}
}

int main(int argc, char *argv[]){
	//parse the parameter
	int c;
	bool sl=false;
	int cutoff=0;
	string fqfilename,outputpath;
	while ((c = getopt(argc,argv,"i:o:c:q")) != -1) {
		switch (c) {
			case 'i':
				fqfilename = optarg;
				break;
			case 'o':
				outputpath = optarg;
				break;
			case 'c':
				cutoff = atoi(optarg);
				break;
			case 'q':
				sl = true;
				break;
			case '?':
				cout << "Unknown option " << optopt << endl;
				break;
			default:
				break;
		}
	}
	if ((!fqfilename.length()) || (!outputpath.length())) {
		cerr << "Input parameters:"<<endl;
		cerr << "pathname:              " << fqfilename << endl;
		cerr << "output file name:    " << outputpath << endl;
		cerr << "sl: " << sl << endl;
		cerr << "cutoff: " << cutoff << endl;
		cerr << "If you are using wildcard, please try quote(\") the string." << endl;
		usage(argv[0]);
		return -1;
	}
	
	//QA filter
	string lowQualitypath="lq"+outputpath;
	string tmppath="tmp"+outputpath;
	if (sl) {
		cout << "Build tmp file" << endl;
		ofstream ofs(tmppath.c_str(),ios_base::out);
		if (!ofs.is_open()) {
			cerr << "Can not open output file." << endl;
			return -1;
		}
			
		ofstream lfs(lowQualitypath.c_str(),ios_base::out);
		if (!lfs.is_open()) {
			cerr << "Can not open output file." << endl;
			return -1;
		}
		
		//read fastq file;
		int lines,badl;
		lines=badl=0;
		map<int,double> pos_qa;
		string s;
		
		gzFile fp;
		kseq_t *seq;
		int l;
		fp = gzopen(fqfilename.c_str(), "r");//open the file handler
		seq = kseq_init(fp);//initialize seq 
		while ((l = kseq_read(seq)) >= 0) {//read sequence  
			lines++;
			if (lines%100000 == 0) {
				cout << '\xd' << "Dealing with reads: " << lines << flush;
			}
			//		cout << "name: " << seq->name.s << endl;
			bool bad=false;
			if (seq->comment.l) {
			//			cout << "comment: " << seq->comment.s << endl;
			stringstream ss(seq->comment.s);
			vector<string> comment;
			while (getline(ss,s,':')) {
					comment.push_back(s);
				}
				if (comment.size()>=3) {//for illumina sequence, the comments is something like 1:Y:18:ATCACG
					if (comment[0]=="2") {//the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
						int paired=1;//no use till now.
					}
					if (comment[1]=="Y") {//Y is the read fails filter (read is bad), N otherwise
						badl++;
						bad=true;
					}
				}
			}		
			string agct(seq->seq.s);
			string qual;
			int max=0;
			int min=127;
			if (seq->qual.l) {
				qual = seq->qual.s;
				if (cutoff>0) {
					getQA(qual,pos_qa,lines,max,min);
					if (max > 74 && min > 58) {
						if ((min - 64) < cutoff) {
							bad=true;
						}
					}else {
						if ((min-33) <  cutoff) {
							bad=true;
						}
					}
				}
			}else {
				qual = "";
			}
			transform(agct.begin(),agct.end(),agct.begin(),::toupper);//convert to upper case
			if(!bad){
				ofs << "@" << seq->name.s;
				if (seq->comment.l) {
					ofs << " " << seq->comment.s;
				}
				ofs << endl << agct << endl;
				if (seq->qual.l) ofs << "+" << endl << qual << endl;
			}else {
				lfs << "@" << seq->name.s;
				if (seq->comment.l) {
					lfs << " " << seq->comment.s;
				}
				lfs << endl << agct << endl;
				if (seq->qual.l) lfs << "+" << endl << qual << endl;			
			}
		}
		kseq_destroy(seq);//destroy seq
		gzclose(fp);//close the file handler
		
		ofs.close();
		lfs.close();	
	}
	//write to s2c files
	cout << "write to s2c files" << endl;
	map<string,int> codemap = initCodeMap();
	map<int,string> rcodemap = mapflip(codemap);
	if (sl) {
		fqfilename=tmppath;
	}
	fq2sc(fqfilename,outputpath,codemap,rcodemap);
	if (sl){
		remove(fqfilename.c_str());
	}
	return 0;
}