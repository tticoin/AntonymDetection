/*
 * WordCounter.cpp
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */

#include "Utils.h"
#include "WordCounter.h"
#include <boost/algorithm/string.hpp>
#include <sstream>
#include <list>
#include <gzstream.h>
#include <map>
#include <cstdio>
using namespace std;

void coin::GoogleNgramWordCounter::count(const string& gzfile){
	igzstream in(gzfile.c_str());
	if(!in){
		cerr << gzfile << " not found!!" << endl;
		return;
	}
	string line;
	while(in && getline(in, line)){
		preprocess(line);
		vector<string> words;
		boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);
		int size = words.size() - 3;
		if(size != 5){
			// TODO: utf-8?
			cerr << "skip line:<" << line << ">" << endl;
			continue;
		}
		long count = atol(words[size+1].c_str());
		for(int i = 0;i < size;i++){
			unsigned int found = words[i].find_first_of("_", 0);
			if(found != string::npos && found != 0){
				words[i] = words[i].substr(0, found);
			}
			if(dict_counts_.find(words[i]) != dict_counts_.end()){
				dict_counts_[words[i]] += count / 5.;
			}else{
				saving_->add(words[i], count / 5.);
			}
		}
	}
	in.close();
}

void coin::WordCounter::set_dict(const string& dict){
	ifstream in(dict.c_str());	
	if(!in){
		cerr << dict << " not found!!" << endl;
		exit(0);
	}
	string line;
	string word;
	while(in && getline(in, line)){
		preprocess(line);
		istringstream iss(line);
		while(iss >> word){
			dict_counts_[word] = 0.;
		}
	}
}


void coin::WordCounter::set_dict(){
	set_dict("dict/synonym.txt");
	set_dict("dict/antonym.txt");
}

void coin::WordCounter::count(const string& gzfile){
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
		while(iss >> word){
			if(dict_counts_.find(word) != dict_counts_.end()){
				dict_counts_[word]++;
			}else{
				saving_->add(word);
			}
		}
	}
	in.close();
}

void coin::WordCounter::write_frequencies(const string& file, int minimum) const{
	ofstream out(file, ios::out|ios::binary|ios::trunc);
	double total = saving_->total();
	out.write((char *)&total, sizeof(double));
    Bucket<string> *bucket = saving_->least_bucket();
    int64_t counter = 0;
    do{
        bucket = bucket->prev();
		if(bucket->count() < minimum){
			bucket = bucket->next();
			break;
		}
        Element<string> *child = bucket->child();
        map<string, Element<string>*> children;
        do{
        	children[child->item()] = child;
            child = child->next();
            ++counter;
        }while(child != bucket->child());
        int childs = children.size();
        out.write((char *)&childs, sizeof(int));
        double count = bucket->count();
        out.write((char *)&count, sizeof(double));
        for(pair<string, Element<string>*> child:children){
        	short size = child.first.size();
         	out.write((char *)&size, sizeof(short));
        	out.write(child.first.c_str(), sizeof(char) * size);
        }
        if(counter > saving_->max_count()){
        	break;
        }
    }while(bucket != saving_->least_bucket());
	for(pair<string, double> dict_count:dict_counts_){
		int childs = 1;
        out.write((char *)&childs, sizeof(int));
        out.write((char *)&dict_count.second, sizeof(double));
		short size = dict_count.first.size();
		out.write((char *)&size, sizeof(short));
		out.write(dict_count.first.c_str(), sizeof(char) * size);
	}
	out.close();
    cerr << counter << " words with at least " << bucket->count() << " frequencies are written" << endl;
}

template <class T>
coin::SpaceSaving<T>::~SpaceSaving(){
	Bucket<T> *bucket = least_bucket_;
	while(bucket->next() != least_bucket_){
		Element<T> *child = bucket->child();
		while(child->next() != bucket->child()){
			Element<T> *prev_child = child;
			child = child->next();
			delete prev_child;
		}
		Bucket<T> *prev_bucket = bucket;
		bucket = bucket->next();
		delete prev_bucket;
	}
};

template <class T>
void coin::SpaceSaving<T>::add(const T& item, double k){
	assert(k > 1e-10);
	total_ += k;
	Element<T> *elem;
	if(least_bucket_ == nullptr){
		++nitems_;
		least_bucket_ = new Bucket<T>(k);
		least_bucket_->set_prev(least_bucket_);
		least_bucket_->set_next(least_bucket_);
		elem = new Element<T>(item);
		least_bucket_->insert(elem);
		counter_[item] = elem;
	}else if(counter_.find(item) == counter_.end()){
		++nitems_;
		assert(nitems_ > 0);
		if(nitems_ > max_buckets_){
			nitems_ = max_buckets_;
			elem = least_bucket_->child();
			double new_count = least_bucket_->count() + k;
			assert(new_count > 0);
			least_bucket_->detach(elem);
			counter_.erase(elem->item());
			elem->set_item(item, least_bucket_->count());
			Bucket<T>* target_bucket = least_bucket_;
			while(true){
				target_bucket = target_bucket->next();
				if(equal(target_bucket->count(), new_count)){
					target_bucket->insert(elem);
					break;
				}else if(target_bucket->count() > new_count || target_bucket == least_bucket_){
					// create new bucket
					Bucket<T>* new_bucket = new Bucket<T>(new_count);
					new_bucket->insert(elem);
					new_bucket->insert_before(target_bucket);
					break;
				}
			}
			if(least_bucket_->empty()){
				Bucket<T> *empty_bucket = least_bucket_;
				least_bucket_ = empty_bucket->next();
				empty_bucket->prev()->set_next(empty_bucket->next());
				empty_bucket->next()->set_prev(empty_bucket->prev());
				delete empty_bucket;
			}
		}else{
			elem = new Element<T>(item);
			double new_count = k;
			Bucket<T>* target_bucket = least_bucket_;
			if(equal(target_bucket->count(), new_count)){
				target_bucket->insert(elem);
			}else if(target_bucket->count() > new_count){
				// create new bucket
				Bucket<T>* new_bucket = new Bucket<T>(new_count);
				new_bucket->insert(elem);
				new_bucket->insert_before(target_bucket);
				least_bucket_ = new_bucket;
			}else{
				while(true){
					target_bucket = target_bucket->next();
					if(equal(target_bucket->count(), new_count)){
						target_bucket->insert(elem);
						break;
					}else if(target_bucket->count() > new_count || target_bucket == least_bucket_){
						// create new bucket
						Bucket<T> *new_bucket = new Bucket<T>(new_count);
						new_bucket->insert(elem);
						new_bucket->insert_before(target_bucket);
						break;
					}
				}
			}
		}
		counter_[item] = elem;
	}else{
		elem = counter_[item];
		Bucket<T> *last_bucket = elem->parent();
		double new_count = last_bucket->count() + k;
		assert(new_count > 0);
		last_bucket->detach(elem);
		Bucket<T>* target_bucket = last_bucket;
		while(true){
			target_bucket = target_bucket->next();
			if(equal(target_bucket->count(), new_count)){
				target_bucket->insert(elem);
				break;
			}else if(target_bucket->count() > new_count || target_bucket == least_bucket_){
				// create new bucket
				Bucket<T> *new_bucket = new Bucket<T>(new_count);
				new_bucket->insert(elem);
				new_bucket->insert_before(target_bucket);
				break;
			}
		}
		if(last_bucket->empty()){
			Bucket<T> *empty_bucket = last_bucket;
			if(last_bucket == least_bucket_){
				least_bucket_ = empty_bucket->next();
			}
			empty_bucket->prev()->set_next(empty_bucket->next());
			empty_bucket->next()->set_prev(empty_bucket->prev());
			delete empty_bucket;
		}
	}
	assert(least_bucket()->prev()->count() >= least_bucket()->count());
	assert(elem->parent());
	assert(elem->item() == item);
	assert(counter_[item] == elem);
}

template <class T>
void coin::Bucket<T>::insert(Element<T> *elem){
    assert(this != nullptr);
    elem->set_parent(this);
    if(child_ == nullptr){
		child_ = elem;
		elem->set_next(elem);
		elem->set_prev(elem);
	}else{
		elem->insert_before(child_);
	}
}

template <class T>
void coin::Bucket<T>::detach(Element<T> *elem){
	if(elem->prev() == elem && elem->next() == elem){
		child_ = nullptr;
	}else{
		if(child_ == elem){
			child_ = elem->next();
		}
		assert(elem->next() != elem);
		elem->prev()->set_next(elem->next());
		assert(elem->prev() != elem);
		elem->next()->set_prev(elem->prev());
		// clean up ...
		elem->set_prev(nullptr);
		elem->set_next(nullptr);
		elem->set_parent(nullptr);
	}
}

template <class I>
void coin::linked_list_item<I>::insert_before(I *next){
	next_ = next;
	prev_ = next->prev();
	prev_->set_next(static_cast<I*>(this));
	next_->set_prev(static_cast<I*>(this));
}


// explicit instantiation of a template
namespace coin{
template class SpaceSaving<string>;
};


