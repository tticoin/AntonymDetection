/*
 * TrainAntonym.cpp
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */

#include "GloveVLBL.h"
#include "Utils.h"
#include <string>
#include <iostream>
#include <cstdlib>
using namespace std;

namespace coin{

template<class T>
void run(const string& freq_file, const string& coocc_file, const string& model_file, int dim, int niter){
	Timer t1;
	GloveVLBL<T> glove(freq_file, coocc_file, dim);
	cerr << "Initialization finished in " << t1.seconds() << " seconds." << endl;

	Timer t2;
	for(int iter = 0;iter < niter;++iter){
		glove.iterate();
	}
	cerr << "Learning finished in " << t2.seconds() << " seconds." << endl;

	Timer t3;
	glove.save_model(model_file);
	cerr << "Model saved to " << model_file << " in " << t3.seconds() << " seconds." << endl;
}


template<class T>
void run(const string& freq_file, const string& coocc_file, const string& model_file, int niter){
    Timer t1;
	GloveVLBL<T> glove(freq_file, coocc_file, model_file);
    cerr << "Initialization finished in " << t1.seconds() << " seconds." << endl;

	Timer t2;
	for(int iter = 0;iter < niter;++iter){
		glove.iterate();
	}
	cerr << "Learning finished in " << t2.seconds() << " seconds." << endl;

	Timer t3;
	glove.save_model(model_file);
	cerr << "Model saved to " << model_file << " in " << t3.seconds() << " seconds." << endl;
}


};

int main(int argc, char *argv[]){
	if(argc < 6){
		cerr << "Usage: " << argv[0] << " freq_file coocc_file model_file dim|- iter" << endl;
		return 0;
	}
	coin::Timer t;
	string freq_file = string(argv[1]);
	string coocc_file = string(argv[2]);
	string model_file = string(argv[3]);
	if(string(argv[4]) == "-"){
		int niter = atoi(argv[5]);
    coin::run<coin::SymmetricWord>(freq_file, coocc_file, model_file, niter);
	}else{
		int dim = atoi(argv[4]);
		int niter = atoi(argv[5]);
    coin::run<coin::SymmetricWord>(freq_file, coocc_file, model_file, dim, niter);
	}
	cerr << "elapsed time: " << t.seconds() << " seconds." << endl;
	return 0;
}
