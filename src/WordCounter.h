/*
 * WordCounter.h
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */
#ifndef WORDCOUNTER_H_
#define WORDCOUNTER_H_

#include <queue>
#include <string>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

namespace coin{

// linked list
template <class I> class linked_list_item{
private:
	I* prev_;
	I* next_;
public:
	linked_list_item():prev_(nullptr), next_(nullptr){}
	virtual ~linked_list_item(){}
	// insert this item before next
	// (next should be in a linked list)
	void insert_before(I* next);
	I* next(){
		return next_;
	}
	void set_next(I* next) {
		next_ = next;
	}
	I* prev(){
		return prev_;
	}
	void set_prev(I* prev) {
	    prev_ = prev;
	}
};


template <class T> class Bucket;
template <class T>
class Element: public linked_list_item<Element<T> > {
friend class Bucket<T>;
private:
	T item_;
	Bucket<T>* parent_;
	double epsilon_;
	void set_parent(Bucket<T>* parent) {
	    parent_ = parent;
	}
public:
	Element(const T &item):item_(item), parent_(nullptr), epsilon_(0){}
	Bucket<T>* parent() const {
		return parent_;
	}
	const T& item() const {
		return item_;
	}
	double epsilon() const{
		return epsilon_;
	}
	// replace word with epsilon-error
	void set_item(const T& item, double epsilon) {
		item_ = item;
		epsilon_ = epsilon;
	}
};

template <class T>
class Bucket: public linked_list_item<Bucket<T> > {
private:
	double count_;
	Element<T>* child_;
public:
 Bucket(double count):count_(count), child_(nullptr){}
	Element<T>* child(){
		return child_;
	}
	double count() const{
		return count_;
	}
	bool empty(){
		return child_ == nullptr;
	}
	// insert elem to this bucket
	void insert(Element<T> *elem);
	// remove elem from this bucket
	void detach(Element<T> *elem);
};
 
// Metwally, Ahmed, Divyakant Agrawal, and Amr El Abbadi. "Efficient computation of frequent and top-k elements in data streams." Database Theory-ICDT 2005. Springer Berlin Heidelberg, 2005. 398-412.
template <class T>
class SpaceSaving final{
private:
	std::unordered_map<T, Element<T>*> counter_;
	int64_t nitems_;
	int64_t max_count_; // # of (empirically) reliable buckets
	int max_buckets_; // # of buckets
	Bucket<T> *least_bucket_; // least hit items
	bool equal(double d1, double d2) const{
		return std::abs(d1 - d2) < 1e-12;
	}
	double total_;
public:
	SpaceSaving(int64_t max_count):
	nitems_(0), max_count_(max_count), max_buckets_(max_count * 10), least_bucket_(nullptr), total_(0.) {}
	~SpaceSaving();
	// add item to stream-summary
	void add(const T& item, double k = 1.);
	// least hit items (prev=>largest hit item)
	Bucket<T>* &least_bucket() {
		return least_bucket_;
	}
	int64_t max_count() const {
		return max_count_;
	}
	double total() const{
		return total_;
	}
};

template <class T>
class Counter{
protected:
	SpaceSaving<T> *saving_;
public:
	virtual ~Counter(){
		if(saving_){
			delete saving_;
		}
	}
	// store counts to space saving
	virtual void count(const std::string& gzfile) = 0;
	// write out frequence words
	virtual void write_frequencies(const std::string &file, int minimum) const = 0;
	// add dictionary words
	virtual void set_dict() = 0;
};

// count frequent words
class WordCounter: public Counter<std::string>{
protected:
	std::unordered_map<std::string, double> dict_counts_;
	using Counter<std::string>::saving_;
	void set_dict(const std::string& dict);
public:
	virtual void count(const std::string& gzfile);
	WordCounter(int max_count = 50000){
		saving_ = new SpaceSaving<std::string>(max_count);
	}
	void write_frequencies(const std::string &file, int minimum) const;
	void set_dict();
};

// count frequent words
class GoogleNgramWordCounter final: public WordCounter{
public:
	void count(const std::string& gzfile);
    GoogleNgramWordCounter(int max_count = 50000):WordCounter(max_count){}
};



} /* namespace coin */


#endif /* WORDCOUNTER_H_ */
