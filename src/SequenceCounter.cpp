/*
 * SequenceCounter.cpp
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */
#include "Utils.h"
#include "SequenceCounter.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <gzstream.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <deque>
using namespace std;
using namespace boost::filesystem;

int coin::CooccurrenceCounter::calc_max_product(double memory_size) const{
	double rlimit = 0.85 * memory_size * 1073741824/(sizeof(double));
	double n = nwords_;
	while(abs(rlimit - n * (log(n) + 0.1544313298)) > 1e-3) n = rlimit / (log(n) + 0.1544313298);
	return (int)n;
}

long long coin::CooccurrenceCounter::calc_excess_limit(double memory_size) const{
	return (long long)(0.15 * memory_size * 1073741824 / (sizeof(std::_Rb_tree_node<std::pair<double, int64_t> >)));
}

std::string coin::CooccurrenceCounter::filename(const std::string &file, int id) const{
	ostringstream os;
	os << file << ".tmp." << id;
	return os.str();
}

void coin::CooccurrenceCounter::count(const string& gzfile, const std::string &outfile){
	igzstream in(gzfile.c_str());
	if(!in){
		cerr << gzfile << " not found!!" << endl;
		return;
	}
	string line;
	string word;
	while(in && getline(in, line)){
		preprocess(line);
		istringstream iss(line);
		vector<int> local_words;
		while(iss >> word){
			if(word_indices_.find(word) != word_indices_.end()){
				local_words.push_back(word_indices_[word]);
			}else{
				local_words.push_back(-1);
			}
		}
		int size = local_words.size();
		// for all pairs
		for(int i = 0;i < size;++i){
			int32_t w1 = local_words[i];
			if(w1 == -1)continue;
			for(int j = i+1;j < size;++j){
				int32_t w2 = local_words[j];
				if(w2 == -1)continue;
				update_count(w1, w2, 1., outfile);
				update_count(w2, w1, 1., outfile);
			}
		}
	}
	in.close();
}

void coin::CooccurrenceCounter::write_excesses(const string& file) const{
	ofstream out(file, ios::out | ios::binary | ios::trunc);
	if(!out){
		cerr << "can not write to " << file << "!!" << endl;
		return;
	}
	out.write((char *)&excesses_, sizeof(int));
	for(pair<int64_t, double> excess_count:excess_counts_){
		out.write((char *)&excess_count.first, sizeof(int64_t));
		out.write((char *)&excess_count.second, sizeof(double));
	}
}

void coin::CooccurrenceCounter::update_count(int32_t w1, int32_t w2, double count, const std::string &outfile){
	total_ += count;
	if(w2 < count_sizes_[w1]){
		counts_[w1][w2] += count;
	}else{
		int64_t l12 = to_scalar(w1, w2);
		if(excess_counts_.find(l12) == excess_counts_.end()){
			excess_counts_[l12] = count;
			excesses_++;
			if(excesses_ > excess_limit_){
				write_excesses(filename(outfile,nexcess_files_));
				nexcess_files_++;
				excesses_ = 0;
				excess_counts_.clear();
			}
		}else{
			excess_counts_[l12] += count;
		}
	}
}

void coin::CooccurrenceCounter::write_cooccurrences(const std::string &file){
	long long cooccurrences = 0;
	ofstream out(file, ios::out | ios::binary | ios::trunc);
	for(int w1 = 0;w1 < nwords_;++w1){
		for(int w2 = 0; w2 < count_sizes_[w1];++w2){
			if(counts_[w1][w2] == 0.)continue;
			out.write((char *)&w1, sizeof(int32_t));
			out.write((char *)&w2, sizeof(int32_t));
			out.write((char *)&counts_[w1][w2], sizeof(double));
			if(cooccurrences % 1000000 == 0){
				cerr << cooccurrences << " pair: "<< word_strings_[w1] << "<=>" << word_strings_[w2] << ":" << counts_[w1][w2] << endl;
			}
			cooccurrences++;
		}
	}
	cerr << cooccurrences << " co-occurrence pairs stored in memory" << endl;
	if(excesses_ > 0){
		write_excesses(filename(file,nexcess_files_));
		nexcess_files_++;
		excesses_ = 0;
		excess_counts_.clear();
	}
	map<int64_t, double> counts;
	ifstream *excesses[nexcess_files_];
	string files[nexcess_files_];
	int sizes[nexcess_files_];
	int64_t current[nexcess_files_];
	for(int i = 0;i < nexcess_files_;++i){
		files[i] = filename(file, i);
		excesses[i] = new ifstream(files[i], ios::in | ios::binary);
		excesses[i]->read((char *)&sizes[i], sizeof(int));
		double count;
		excesses[i]->read((char *)&current[i], sizeof(int64_t));
		excesses[i]->read((char *)&count, sizeof(double));
		if(counts.find(current[i]) != counts.end()){
			counts[current[i]] += count;
		}else{
			counts[current[i]] = count;
		}
		sizes[i]--;
	}
	while(nexcess_files_ > 0){
		int64_t index = counts.begin()->first;
		double value = counts.begin()->second;
		counts.erase(index);
		pair<int32_t, int32_t> w = to_pair(index);
		out.write((char *)&w.first, sizeof(int32_t));
		out.write((char *)&w.second, sizeof(int32_t));
		out.write((char *)&value, sizeof(double));
		if(cooccurrences % 1000000 == 0){
			cerr << cooccurrences << " pair: "<< word_strings_[w.first] << "<=>" << word_strings_[w.second] << ":" << value << endl;
		}
		cooccurrences++;
		for(int i = 0;i < nexcess_files_;++i){
			if(index == current[i] && sizes[i] != 0){
				double count;
				excesses[i]->read((char *)&current[i], sizeof(int64_t));
				excesses[i]->read((char *)&count, sizeof(double));
				if(counts.find(current[i]) != counts.end()){
					counts[current[i]] += count;
				}else{
					counts[current[i]] = count;
				}
				sizes[i]--;
			}
		}
		bool need_check = false;
		do{
			need_check = false;
			for(int i = 0;i < nexcess_files_;++i){
				if(sizes[i] == 0){
					excesses[i]->close();
					remove(files[i].c_str());
					delete excesses[i];
					if(nexcess_files_ - 1 != i){
						excesses[i] = excesses[nexcess_files_ - 1];
						files[i] = files[nexcess_files_ - 1];
						sizes[i] = sizes[nexcess_files_-1];
						current[i] = current[nexcess_files_ - 1];
						need_check = true;
					}
					nexcess_files_--;
					break;
				}
			}
		}while(need_check);
	}
	for(pair<int64_t, double> count:counts){
		int64_t index = count.first;
		double value = count.second;
		pair<int32_t, int32_t> w = to_pair(index);
		out.write((char *)&w.first, sizeof(int32_t));
		out.write((char *)&w.second, sizeof(int32_t));
		out.write((char *)&value, sizeof(double));
		if(cooccurrences % 1000000 == 0){
			cerr << cooccurrences << " pair: "<< word_strings_[w.first] << "<=>" << word_strings_[w.second] << ":" << value << endl;
		}
		cooccurrences++;
	}
	cerr << cooccurrences << " pairs in total ";
	out.close();
}

void coin::CooccurrenceCounter::shuffle(const std::string &file) const{
	char *buf_in = new char[BUFSIZE];
	char *buf_out = new char[BUFSIZE];
	string in_filename = file;
	string tmp_filename = in_filename+".tmp";
	ifstream in_file(in_filename, ios::in | ios::binary);
	ofstream out_file(tmp_filename, ios::out | ios::binary);
	in_file.rdbuf()->pubsetbuf(buf_in, BUFSIZE);
	out_file.rdbuf()->pubsetbuf(buf_out, BUFSIZE);
	vector<CooccEntry> entries;
	entries.reserve(SHUFFLE_SIZE+1); 
	int nentries = 0;
	while(true){
		int32_t w1, w2;
		double count;
		in_file.read((char *)&w1, sizeof(int32_t));
		if(in_file.eof())break;
		in_file.read((char *)&w2, sizeof(int32_t));
		in_file.read((char *)&count, sizeof(double));
		entries.push_back(CooccEntry(w1, w2, count));
		nentries++;
		if(nentries > SHUFFLE_SIZE){
			std::shuffle(entries.begin(), entries.end(), default_random_engine(0));
			for(CooccEntry entry:entries){
				out_file.write((char*)&entry.w1, sizeof(int32_t));
				out_file.write((char*)&entry.w2, sizeof(int32_t));
				out_file.write((char*)&entry.count, sizeof(double));
			}
			entries.clear();
			nentries = 0;
		}
	}
	in_file.close();
	if(nentries > 0){
		std::shuffle(entries.begin(), entries.end(), default_random_engine(0));
		for(CooccEntry entry:entries){
			out_file.write((char*)&entry.w1, sizeof(int32_t));
			out_file.write((char*)&entry.w2, sizeof(int32_t));
			out_file.write((char*)&entry.count, sizeof(double));
		}
	}
	out_file.close();
	delete[] buf_in;
	delete[] buf_out;
	try{
		remove(path(in_filename));
		rename(path(tmp_filename), path(in_filename));
	}catch(filesystem_error& ex) {
		cerr << ex.what() << endl;
		exit(-1);
	}
}

void coin::SequenceCounter::count(const string& gzfile, const std::string &outfile){
	igzstream in(gzfile.c_str());
	if(!in){
		cerr << gzfile << " not found!!" << endl;
		return;
	}
	string line;
	string word;
	while(in && getline(in, line)){
		preprocess(line);
		istringstream iss(line);
		vector<int> local_words;
		while(iss >> word){
			if(word_indices_.find(word) != word_indices_.end()){
				local_words.push_back(word_indices_[word]);
			}else{
				local_words.push_back(-1);
			}
		}
		int size = local_words.size();
		for(int i = 0;i < size;i++){
			int32_t w1 = local_words[i];
			if(w1 == -1)continue;
			for(int j = i+1;j < size && j <= i+nwindow_;j++){
				int32_t w2 = local_words[j];
				if(w2 == -1)continue;
				double c = 1.;
				if(weight() == 1){ 
					c = 1. / (j - i);
				}else if(weight() == 2){
					c = 1. + (1. - (j - i)) / nwindow_;
				}
				update_count(w1, w2, c, outfile);
				if(symmetric()){
					update_count(w2, w1, c, outfile);
				}
			}
		}
	}
	in.close();
}
