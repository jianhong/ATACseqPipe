/*
 *  fq2sc.h
 *  
 *
 *  Created by Jianhong Ou on 3/28/12.
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

#ifndef FQ2SC_H_H
#define FQ2SC_H_H

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <tr1/unordered_map>
#include <tr1/memory>
using namespace std;
using namespace std::tr1;

typedef unordered_map<string,int> s2cmap;
typedef pair<string,int> s2cpair;
struct IntCmp {
    bool operator()(const s2cpair &lhs, const s2cpair &rhs) {
        return lhs.second > rhs.second;
    }
};
map<string,int> initCodeMap(){
	map<string,int> codemap;
	int i=0;
	/* map: AGCTNS(stop) use*/
	char coden[]="AGCTNS";
	for (int w=0; w<6; w++) {
		for (int m=0; m<6; m++) {
			for (int n=0; n<6; n++) {
				string s="";
				s+=coden[w];
				s+=coden[m];
				s+=coden[n];
				codemap[s]=i++;
			}
		}
	}
	return codemap;
}

string encode(string s, map<string,int> codemap){
	string x="";
	for (int i=0; i<s.length(); i=i+3) {
		string ss=s.substr(i,3);
		if (ss.length()<3) {
			ss+="SSS";
			ss=ss.substr(0,3);
		}
		x+=(char)(codemap[ss]-128);//char is signed int, can only -128<= char <=127;
	}
	return x;
}
string decode(string s, map<int,string> rcodemap){
	string x="";
	for (int i=0; i<s.length(); i++) {
		x+=rcodemap[(int)(s[i]+128)];
	}
	size_t found=x.find('S');
	if(found!=string::npos){
		x=x.substr(0,found);
	}
	return x;
}
map<int,string> mapflip(map<string,int> codemap){
	map<int,string> rcodemap;
	for (map<string,int>::iterator it=codemap.begin(); it!=codemap.end(); it++) {
		rcodemap[it->second] = it->first;
	}
	return rcodemap;
}

void fq2sc(string infilename, string outfilename, map<string,int> codemap, map<int,string> rcodemap){
//	map<string,int> codemap = initCodeMap();
//	map<int,string> rcodemap = mapflip(codemap);
/*	for (map<string,int>::iterator it=codemap.begin(); it!=codemap.end(); it++) {
		cout << it->first << "\t" << it->second << endl;
	}
	for (map<int,string>::iterator it=rcodemap.begin(); it!=rcodemap.end(); it++) {
		cout << it->first << "\t" << it->second << endl;
	}*/
	s2cmap s2c;
	ifstream ifs(infilename.c_str(), ios_base::in);
//	ofstream log("log.txt", ios_base::out);
	string s;
	bool at=false;
	bool plus=false;
	int lines=0;
	if (ifs.is_open()) {
		while (!ifs.eof()) {
			if(getline(ifs,s)){
				lines++;
//				if (lines%1000 == 0) {
//					cout << '\xd' << "Dealing with lines: " << lines << flush;
//				}
				if(at){
//					log << s << "\t" << atoi(encode(s,codemap).c_str()) << "\t" << decode(encode(s,codemap),rcodemap) << endl << flush;
					s2c[encode(s,codemap)]++;
					at=false;
				}else if(s[0]=='@' && (!plus)){
					at=true;
					continue;
				}else if (s[0]=='+' && (!plus)) {
					plus=true;
					continue;
				}else if (plus) {
					plus=false;
					continue;
				}else {
					cerr << "err at : " << s << "\t" << lines << endl;
				}

			}
		}
//		cout << endl;
		ifs.close();
		//sort the hash
		vector<s2cpair> s2cvec(s2c.begin(),s2c.end());
		sort(s2cvec.begin(),s2cvec.end(),IntCmp());
		//output the hash
		ofstream ofs(outfilename.c_str(),ios_base::out);
		if (ofs.is_open()) {
			for (vector<s2cpair>::iterator it=s2cvec.begin(); it!=s2cvec.end(); it++) {
				ofs << decode(it->first,rcodemap) << "\t" << it->second << endl;
			}
			ofs.close();
		}else {
			cerr << "Can not open " << outfilename << endl;
		}
	}else {
		cerr << "Can not open " << infilename << endl;
	}
//	log.close();
}

#endif

/*	this maybe used later, this time we do not use it.
 *	first cut as 3 char per unit
 *	make first char record the length (max length:256)
 *	from the second char to record the sequence
 *	every char is 1byte, 8bit. will recode 3 base
 *	A: 00	G: 01	C: 10	T: 11
 *	the Highest 2bit record where is the first N
 *	H: 8 7 6 5 4 3 2 1 :L
 *	8 7: 00, no N	01, N in 6 5	10, N in 4 3	11, N in 2 1
 *	If there is N in 6 5 or 4 3, in 6 5 or 4 3 will record the next N
 *	6 5: 00, no N	01, N in 4 3	10, N in 2 1
 *	4 3: 00, no N	01, N in 2 1
 *	map sample:
 *		AAA	00 00 00 00
 *		AAG	00 00 00 01
 *		AAN 11 00 00 00
 *		ANA	10 00 00 00
 *		NAA	01 00 00 00
 *		ANN	10 00 01 00
 *		NAN	01 10 00 00
 *		NNA	01 01 00 00
 *		NNN	01 01 01 00
 */