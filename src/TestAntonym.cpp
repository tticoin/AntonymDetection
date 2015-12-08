/*
 * GloveVLBLAntonym.cpp
 *
 *  Created on: 2014/09/18
 *      Author: miwa
 */

#include "GloveVLBL.h"
#include <string>
#include <iostream>
using namespace std;
using namespace coin;

template<class T>
void run(string& freq_file, string& model_file, string& test_set){
	GloveVLBL<T> glove(freq_file, model_file);
	cerr << "Model loaded." << endl;
	glove.antonym_score(test_set);
}

int main(int argc, char *argv[]){
	if(argc < 4){
		cerr << "Usage: " << argv[0] << " freq_file model_file test_set" << endl;
		return 0;
	}
	string freq_file = string(argv[1]);
	string model_file = string(argv[2]);
	string test_set = string(argv[3]);
  run<SymmetricWord>(freq_file, model_file, test_set);
	return 0;
}
