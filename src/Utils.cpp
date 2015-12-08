/*
 * Utils.cpp
 *
 *  Created on: 2014/09/14
 *      Author: miwa
 */
#include "Utils.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cassert>
using namespace std;

int normalize_digit(int c){
	if(isdigit(c)){
		return '7';
	}
	return c;
}

void coin::preprocess(string& str){
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	//transform(str.begin(), str.end(), str.begin(), normalize_digit);
}

double coin::load_word_counts(const string& file, unordered_map<string, int32_t> &word_indices, vector<string> &word_strings, vector<double> &word_counts){
	ifstream in(file, ios::in | ios::binary);
	if(!in){
		cerr << file << " not found!!" << endl;
		exit(0);
	}
	int index = 0;
	double all_total;
	// ignore total counts...
	in.read((char *)&all_total, sizeof(double));
	double total = 0.;
	while(!in.eof()){
		int childs;
		double count;
		in.read((char *)&childs, sizeof(int));
		if(in.eof())break;
		in.read((char *)&count, sizeof(double));
		total += childs * count;
		while(childs-- > 0){
			short size;
			in.read((char *)&size, sizeof(short));
			char s[size+1];
			in.read((char *)&s[0], sizeof(char) * size);
			s[size] = '\0';
			string word(s);
			assert(word_indices.find(word) == word_indices.end());
			word_indices[word] = index;
			word_strings.push_back(word);
			word_counts.push_back(count);
			index++;
		}
	}
	in.close();
	//cerr << index << " words starting from \"" << word_strings[0] << "\" read" << endl;
	return total;
}

double coin::cosine(const Vector& l, const Vector& r){
	return l.dot(r) / l.norm() / r.norm();
}

pair<int32_t, int32_t> coin::to_pair(int64_t i){
	return pair<int32_t, int32_t>(i >> 32, i & (((int64_t)2 << 32) - 1));
}

int64_t coin::to_scalar(int32_t i1, int32_t i2){
	return (((int64_t)i1) << 32) | i2;
}

double coin::sigmoid(const double x){
	return 1. / (1. + exp(-x));
}

void coin::print_mem(){
	ifstream proc_stream("/proc/self/statm");
	long long VmSize = 0, VmRSS = 0, Share = 0;
	proc_stream >> VmSize >> VmRSS >> Share;
	proc_stream.close();
	cerr << "Memory usage: " << VmRSS  * getpagesize () / (1024.0 * 1024.0) << " MiB." << endl;
}
