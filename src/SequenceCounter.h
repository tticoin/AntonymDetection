/*
 * SequenceCounter.h
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */

#ifndef SEQUENCECOUNTER_H_
#define SEQUENCECOUNTER_H_

#include "Utils.h"
#include <cassert>
#include <unordered_map>
#include <map>
#include <string>
#include <vector>

namespace coin {

class CooccurrenceCounter {
protected:
	static constexpr int BUFSIZE = 1024*1024; 
	static constexpr int SHUFFLE_SIZE = 1000000;
	// members for frequent words
	std::unordered_map<std::string, int32_t> word_indices_; // string => index
	std::vector<std::string> word_strings_; // index => string
	std::vector<double> word_counts_; // index => count
	int nwords_; // # of words
	double total_; // # of word occurrences

	// members for storing co-occurrences in memory
	int max_product_;
	int *count_sizes_;
	double **counts_;
	int calc_max_product(double memory_size) const;

	// members for storing co-occurrences in files
	int nexcess_files_ = 0; // # of flushed files
	long long excess_limit_; // limit to flush
	long long excesses_ = 0; // counter for excess counts (excesses_ == excess_counts_.size())
	std::map<int64_t, double> excess_counts_; // exceeded co-occurrences
	long long calc_excess_limit(double memory_size) const;
	void write_excesses(const std::string& file) const;
	std::string filename(const std::string &file, int id) const;

	void update_count(int32_t w1, int32_t w2, double count, const std::string &outfile);
public:
	CooccurrenceCounter(const std::string& word_counts, double memory_size = 10.): total_(0.){// 100MB
		load_word_counts(word_counts, word_indices_, word_strings_, word_counts_);
		nwords_ = word_indices_.size();
		count_sizes_ = new int[nwords_];
		counts_ = new double*[nwords_];
		max_product_ = calc_max_product(memory_size);
		for(int n = 0;n < nwords_;n++){
			count_sizes_[n] = max_product_ / (n+1);
			counts_[n] = new double[count_sizes_[n]]();
			assert(count_sizes_[n] == 0 || counts_[n][0] == 0.);
		}
		excess_limit_ = calc_excess_limit(memory_size);
	}
	virtual ~CooccurrenceCounter(){
		for(int n = 0;n < nwords_;++n){
			delete[] counts_[n];
		}
		delete[] counts_;
		delete[] count_sizes_;
	}
	virtual void count(const std::string& gzfile, const std::string &outfile);
	// write out co-occurred pairs
	virtual void write_cooccurrences(const std::string &file);
	void shuffle(const std::string &file) const;
	double total() const{return total_;}
};

class SequenceCounter : public CooccurrenceCounter {
protected:
	int nwindow_;
	bool symmetric_;
	int weight_;
public:
 SequenceCounter(const std::string& word_counts, int nwindow = 10, double memory_size = 10., bool symmetric = false, int weight = 1):
   CooccurrenceCounter(word_counts, memory_size), nwindow_(nwindow), symmetric_(symmetric), weight_(weight){}
	virtual void count(const std::string& gzfile, const std::string &outfile);
	bool symmetric() const{return symmetric_;}
	int weight() const{return weight_;}
};


class CooccEntry {
public:
	int32_t w1;
	int32_t w2;
	double count;
 CooccEntry(int32_t w1_, int32_t w2_, double count_):w1(w1_), w2(w2_), count(count_){}
};
 

} /* namespace coin */

#endif /* SEQUENCECOUNTER_H_ */
