/*
 * SequenceCount.cpp
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */

#include "SequenceCounter.h"
#include <iostream>
#include <boost/filesystem.hpp>
using namespace coin;
using namespace std;
using namespace boost::filesystem;

// count co-occurrences in windows
int main(int argc, char *argv[]){
	if(argc < 4){
		cerr << "Usage: " << argv[0] << " dir freq_file seq_file window weight" << endl;
		cerr << "weight\t2: ivlbl, 1: glove, 0: noweight (default: 1)" << endl;
		return 0;
	}
	int window = atoi(argv[4]); // check [-window, window]
	int weight = 1;
	if(argc > 5){
		weight = atoi(argv[5]);
	}
	cerr << "weight: " << weight << endl; 
	SequenceCounter counter(argv[2], window, 10, 1, weight);
	path p(argv[1]);
	if(exists(p) && is_directory(p)){
	    for(directory_iterator it(p); it != directory_iterator(); ++it){
	    	if(is_regular_file(*it) && (extension(*it) == ".gz")){
	    		cerr << "reading " << it->path().string() << endl;
	    		counter.count(it->path().string(), argv[3]);
	    	}
	    }
	}
	counter.write_cooccurrences(argv[3]);
	counter.shuffle(argv[3]);
	return 0;
}
