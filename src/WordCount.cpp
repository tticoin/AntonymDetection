/*
 * WordCount.cpp
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */

#include "WordCounter.h"
#include <iostream>
#include <boost/filesystem.hpp>
using namespace coin;
using namespace std;
using namespace boost::filesystem;

int main(int argc, char *argv[]){
	if(argc < 4){
		cerr << "Usage: " << argv[0] << " dir count file [minimum_freq] [use_dict]" << endl;
		cerr << "use_dict: 1 or 0" << endl;
		return 0;
	}
	int count = atoi(argv[2]);
	assert(count > 0);
	WordCounter counter(count);
	path p(argv[1]);
	int use_dictionary = 0;
	if(argc > 5){
		use_dictionary = atoi(argv[5]);
		if(use_dictionary){
			counter.set_dict();
		}
	}
	if(exists(p) && is_directory(p)){
	    for(directory_iterator it(p); it != directory_iterator(); ++it){
	    	if(is_regular_file(*it) && (extension(*it) == ".gz")){
	    		cerr << "reading " << it->path().string() << endl;
	    		counter.count(it->path().string());
	    	}
	    }
	}
	
	int minimum = 0;
	if(argc > 4){
		minimum = atoi(argv[4]);
	}
	counter.write_frequencies(argv[3], minimum);
	return 0;
}

