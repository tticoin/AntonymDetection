/*
 * GloveVLBL.cpp
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */
#include "GloveVLBL.h"
#include <gzstream.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <unordered_set>
#include <cmath>
using namespace std;
using namespace boost::filesystem;

template<class T>
coin::GloveVLBL<T>::GloveVLBL(const string& freq_file, const string& coocc_file, int dim): coocc_file_(coocc_file), dim_(dim) {
	total_ = load_word_counts(freq_file, word_indices_, word_strings_, word_counts_);
	function<double()> rand(bind(uniform_real_distribution<double>(-.5/dim, .5/dim), default_random_engine(0)));
	words_.resize(word_indices_.size());
	int nwords = 0;
	for(pair<string, int> word:word_indices_){
		words_[word.second] = new T(word.second, dim_, rand);
		nwords++;
	}
	assert(nwords == word_indices_.size());
	init_table();
}

template<class T>
coin::GloveVLBL<T>::GloveVLBL(const string& freq_file, const string& model_file) : coocc_file_("") {
	if(freq_file == "-g"){
		load_glove_vector(model_file, false);
	}else if(freq_file == "-l"){
		// log-bilinear
		load_glove_vector(model_file, true);
	}else{
		total_ = load_word_counts(freq_file, word_indices_, word_strings_, word_counts_);
		load_model(model_file);
	}
}

template<class T>
coin::GloveVLBL<T>::GloveVLBL(const string& freq_file, const string& coocc_file, const string& model_file): coocc_file_(coocc_file) {
	total_ = load_word_counts(freq_file, word_indices_, word_strings_, word_counts_);
	load_model(model_file);
	init_table();
}

template<class T>
bool more(const pair<T, double>& left, const pair<T, double>& right ) {
    if( left.second < right.second ){
    	return false;
    }else if(right.second < left.second){
    	return true;
    }else{
    	return left.first < right.first;
    }
}



template<class T>
string coin::GloveVLBL<T>::filename(int index) const{
	ostringstream oss;
	oss << tmpdir_ << "/" << index;
	return oss.str();
}

string current(){
	ostringstream oss;
	chrono::system_clock::time_point p = chrono::system_clock::now();
	time_t t = chrono::system_clock::to_time_t(p);
  char time_string[100];
  struct tm *tm = localtime(&t);
  strftime(time_string,sizeof(time_string),"%F_%T",tm);
  oss << time_string;
	return oss.str();
}

template<class T>
void coin::GloveVLBL<T>::load_dict(const std::string& dict, std::unordered_map<int, set<int>>& map){
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
		iss >> word;
		if(word_indices_.find(word) != word_indices_.end()){
			int target_id = word_indices_[word];
			while(iss >> word){
				if(word_indices_.find(word) != word_indices_.end()){
					int ref_id = word_indices_[word];
					if(target_id == ref_id)continue;
					if(target_id < ref_id){
						if(map.find(target_id) == map.end()){
							map[target_id] = set<int>();
						}					
						map[target_id].insert(ref_id);
					}else{
						if(map.find(ref_id) == map.end()){
							map[ref_id] = set<int>();
						}					
						map[ref_id].insert(target_id);
					}
				}
			}
		}
	}
}

template<class T>
void coin::GloveVLBL<T>::load_dict(){
	load_dict("data/dict/synonym.txt", synonyms_);
	load_dict("data/dict/antonym.txt", antonyms_);
	cerr << "synonyms:" << synonyms_.size() << ", " << "antonyms:" << antonyms_.size() << endl;
}

template<class T>
bool coin::GloveVLBL<T>::includes_pair(unordered_map<int, set<int> >& map, int w1, int w2){
	if(w1 < w2){
		if(map.find(w1) != map.end() && map[w1].find(w2) != map[w1].end()){
			return true;
		}
	}else{
		if(map.find(w2) != map.end() && map[w2].find(w1) != map[w2].end()){
			return true;
		}
	}
	return false;
}

template<class T>
void coin::GloveVLBL<T>::init_table(){
	load_dict();
	nwords_ = words_.size();
    tmpdir_ = "tmp/"+current();
	if(!is_directory(tmpdir_)){
		if(!create_directories(tmpdir_)){
			cerr << "cannot create directory: " << tmpdir_ << endl;
			exit(-1);
		}		
	}
	vector<pair<int32_t, count_t>> neg_sample_word_counts;	
	neg_sample_word_counts.reserve(nwords_);
	double neg_sample_total = 0.;
	vector<double> subsampling;
	subsampling.resize(nwords_);
	int index = 0;
	for(double wc:word_counts_){
		//subsampling (word2vec)
		double s = wc / (SUBSAMPLING * total_);
		double p = min((sqrt(s) + 1.) / s, 1.);
		//double p = min(sqrt(total_ * SUBSAMPLING / wc), 1.); // mikolov's paper
		subsampling[index] = p;
		//neg_sampling
		count_t pc = pow(wc, .75);
		neg_sample_word_counts.push_back(pair<int32_t, count_t>(index, pc));
		neg_sample_total += pc;
		index++;
	}
	vector<count_t> subsampled_counts;
	subsampled_counts.resize(nwords_, 0);
	{
		ifstream coocc_in(coocc_file_, ios::in | ios::binary);
		if(!coocc_in){
			cerr << coocc_file_ << " not found!!" << endl;
			exit(0);
		}
		int32_t w1, w2;
		double count;
		double subsampled_count_total = 0.;
		int64_t count_total = 0;
		int64_t total = 0;
		while(!coocc_in.eof()){
			coocc_in.read((char *)&w1, sizeof(int32_t));
			if(coocc_in.eof())break;
			coocc_in.read((char *)&w2, sizeof(int32_t));
			coocc_in.read((char *)&count, sizeof(double));
			
			double subsampled_count = count * subsampling[w1] * subsampling[w2];
			subsampled_counts[w1] += subsampled_count;
			subsampled_count_total += subsampled_count;
			count_total += count;
			total++;
		}
		coocc_in.close();
		cerr << count_total << " (" << (int64_t) subsampled_count_total << " sampled) cooccurrences, " << total << " unique cooccurrences" << endl;
	}
	set<int64_t> antonym_set;
	for(pair<int, std::set<int> > pairs:antonyms_){
		int w1 = pairs.first;
		for(int w2: pairs.second){
			antonym_set.insert(to_scalar(w2, w1));
			antonym_set.insert(to_scalar(w1, w2));
		}
	}
	set<int64_t> synonym_set;
	for(pair<int, std::set<int> > pairs:synonyms_){
		int w1 = pairs.first;
		for(int w2: pairs.second){
			synonym_set.insert(to_scalar(w2, w1));
			synonym_set.insert(to_scalar(w1, w2));
		}
	}
	sort(neg_sample_word_counts.begin(), neg_sample_word_counts.end(), more<int32_t>);
	{
		char *buf[FILES];
		ofstream *tmp_files[FILES];
		for(int i = 0;i < FILES;++i){
			buf[i] = new char[BUFSIZE];
			tmp_files[i] = new ofstream(filename(i), ios::out | ios::binary | ios::trunc);
			tmp_files[i]->rdbuf()->pubsetbuf(buf[i], BUFSIZE);
		}
   		
		unsigned long long next_random = 0.;

		int64_t updates = 0;
		double total_counts = 0., total_neg_counts = 0.;
		int64_t antonym_matches = 0, synonym_matches = 0;
		int file_idx = 0;
		vector<int32_t> *cooccurrences = new vector<int32_t>[nwords_];
		{	
			ifstream coocc_in(coocc_file_, ios::in | ios::binary);
			if(!coocc_in){
				cerr << coocc_file_ << " not found!!" << endl;
				exit(0);
			}
			int32_t w1, w2;
			double coocc_count;
			while(!coocc_in.eof()){
				coocc_in.read((char *)&w1, sizeof(int32_t));
				if(coocc_in.eof())break;
				coocc_in.read((char *)&w2, sizeof(int32_t));
				coocc_in.read((char *)&coocc_count, sizeof(double));

				assert(w1 < nwords_);
				cooccurrences[w1].push_back(w2);

				count_t count = coocc_count * subsampling[w1] * subsampling[w2];
				count_t cnt = (subsampled_counts[w1] * K) / neg_sample_total;
				count_t neg_count = cnt * pow(word_counts_[w2], .75);
				if(count < MIN_COUNT){
					count = 0.;
				}else if(count < 1.){
					next_random = next_random * (unsigned long long)25214903917 + 11;
					count = ((next_random & 0xFFFF) / (count_t)65536 < count) ? 1. : 0.;
				} 
				if(neg_count < MIN_COUNT){
					neg_count = 0.;
				}else if(neg_count < 1.){
					next_random = next_random * (unsigned long long)25214903917 + 11;
					neg_count = ((next_random & 0xFFFF) / (count_t)65536 < neg_count) ? 1. : 0.;
				} 
				if(count > 0. || neg_count > 0.){
					if(synonym_set.find(to_scalar(w1, w2)) != synonym_set.end()){
						continue;
					}
					if(antonym_set.find(to_scalar(w1, w2)) != antonym_set.end()){
						continue;
					}
					total_counts += count;
					total_neg_counts += neg_count;
					tmp_files[file_idx]->write((char *)&w1, sizeof(int32_t));
					tmp_files[file_idx]->write((char *)&w2, sizeof(int32_t));
					tmp_files[file_idx]->write((char *)&count, sizeof(count_t));
					tmp_files[file_idx]->write((char *)&neg_count, sizeof(count_t));
					updates++;
					file_idx = (file_idx+1) % FILES;
				}
			}
			coocc_in.close();
		}
		for(int i = 0;i < nwords_;++i){
			pair<int32_t, count_t>& w1cnt = neg_sample_word_counts[i];

			vector<int32_t>& coocc = cooccurrences[w1cnt.first];
			sort(coocc.begin(), coocc.end());

			int l = 0, r = nwords_ - 1, m = (l+r) / 2;
			count_t cnt = (subsampled_counts[w1cnt.first] * K) / neg_sample_total;
			count_t icnt = MIN_COUNT / cnt;
			while(l < r){
				if(neg_sample_word_counts[m].second >= icnt){
					l = m + 1;
				}else{
					r = m;
				}
				m = (l + r) / 2;
			}
			assert(m == 0 || neg_sample_word_counts[m - 1].second >= icnt);
			assert(m == nwords_ - 1 || neg_sample_word_counts[m].second < icnt);
			for(int j = 0;j < m;++j){
				pair<int32_t, count_t>& w2cnt = neg_sample_word_counts[j];
				if(binary_search(coocc.begin(), coocc.end(), w2cnt.first)){
					continue;
				}
				count_t count = 0.;
				count_t neg_count = cnt * w2cnt.second;

				assert(neg_count >= MIN_COUNT);
				if(neg_count < 1.){
					next_random = next_random * (unsigned long long)25214903917 + 11;
					neg_count = ((next_random & 0xFFFF) / (count_t)65536 < neg_count) ? 1. : 0.;
				} 
				if(neg_count > 0.){
					if(synonym_set.find(to_scalar(w1cnt.first, w2cnt.first)) != synonym_set.end()){
						continue;
					}
					if(antonym_set.find(to_scalar(w1cnt.first, w2cnt.first)) != antonym_set.end()){
						continue;
					}
					total_counts += count;
					total_neg_counts += neg_count;
					tmp_files[file_idx]->write((char *)&w1cnt.first, sizeof(int32_t));
					tmp_files[file_idx]->write((char *)&w2cnt.first, sizeof(int32_t));
					tmp_files[file_idx]->write((char *)&count, sizeof(count_t));
					tmp_files[file_idx]->write((char *)&neg_count, sizeof(count_t));
					updates++;
					file_idx = (file_idx+1) % FILES;
				}
			}
		}
		for(int64_t synonym:synonym_set){
			pair<int32_t, int32_t> p = to_pair(synonym);
			count_t count = SYN_WEIGHT;
			count_t neg_count = 0.;
			synonym_matches++;
			if(antonym_set.find(synonym) != antonym_set.end()){
				continue;
			}
			tmp_files[file_idx]->write((char *)&p.first, sizeof(int32_t));
			tmp_files[file_idx]->write((char *)&p.second, sizeof(int32_t));
			tmp_files[file_idx]->write((char *)&count, sizeof(count_t));
			tmp_files[file_idx]->write((char *)&neg_count, sizeof(count_t));
			file_idx = (file_idx+1) % FILES;
		}
		for(int64_t antonym:antonym_set){
			pair<int32_t, int32_t> p = to_pair(antonym);
			count_t count = 0.;
			count_t neg_count = ANT_WEIGHT;
			antonym_matches++;
			tmp_files[file_idx]->write((char *)&p.first, sizeof(int32_t));
			tmp_files[file_idx]->write((char *)&p.second, sizeof(int32_t));
			tmp_files[file_idx]->write((char *)&count, sizeof(count_t));
			tmp_files[file_idx]->write((char *)&neg_count, sizeof(count_t));
			file_idx = (file_idx+1) % FILES;
		}
		print_mem();
		cerr << "Total updates: " << updates << ", counts: " << total_counts << ", neg_counts: " << total_neg_counts << endl;
		cerr << "Total synonym matches: " << synonym_matches << ", antonym matches: " << antonym_matches << endl;
		for(int i = 0;i < FILES;++i){
			tmp_files[i]->close();
			delete tmp_files[i];
			delete[] buf[i];
		}
		delete[] cooccurrences;
	}
	merge_files();
}


template<class T>
void coin::GloveVLBL<T>::merge_files(){
#ifdef _OPENMP
#pragma omp parallel num_threads(N_THREADS)
#endif
	{
		int size = FILES / N_THREADS;
		int start = omp_get_thread_num() * size;
		char *buf_in = new char[BUFSIZE];
		char *buf_out = new char[BUFSIZE];
		string tmp_filename = filename(omp_get_thread_num())+".tmp";
		ofstream out_file(tmp_filename, ios::out | ios::binary);
		default_random_engine rand(omp_get_thread_num());
		out_file.rdbuf()->pubsetbuf(buf_out, BUFSIZE);
		int nupdates = 0;
		vector<UpdateEntry> updates;
		updates.reserve(SHUFFLE_SIZE+1); 
		for(int i = 0;i < size;++i){
			string in_filename = filename(start+i);
			ifstream in_file(in_filename, ios::in | ios::binary);
			in_file.rdbuf()->pubsetbuf(buf_in, BUFSIZE);
			while(true){
				int32_t w1, w2;
				count_t count, neg_count;
				in_file.read((char *)&w1, sizeof(int32_t));
				if(in_file.eof())break;
				in_file.read((char *)&w2, sizeof(int32_t));
				in_file.read((char *)&count, sizeof(count_t));
				in_file.read((char *)&neg_count, sizeof(count_t));
				updates.push_back(UpdateEntry(w1, w2, count, neg_count));
				nupdates++;
				if(nupdates > SHUFFLE_SIZE){
					shuffle(updates.begin(), updates.end(), rand);
					for(UpdateEntry entry:updates){
						out_file.write((char*)&entry.w1, sizeof(int32_t));
						out_file.write((char*)&entry.w2, sizeof(int32_t));
						out_file.write((char*)&entry.count, sizeof(count_t));
						out_file.write((char*)&entry.neg_count, sizeof(count_t));
					}
					updates.clear();
					nupdates = 0;
				}
			}
			in_file.close();
			try{
				remove(path(in_filename));	
			}catch(filesystem_error& ex) {
				cerr << ex.what() << endl;
				exit(-1);
			}
		}
		if(nupdates > 0){
			shuffle(updates.begin(), updates.end(), rand);
			for(UpdateEntry entry:updates){
				out_file.write((char*)&entry.w1, sizeof(int32_t));
				out_file.write((char*)&entry.w2, sizeof(int32_t));
				out_file.write((char*)&entry.count, sizeof(count_t));
				out_file.write((char*)&entry.neg_count, sizeof(count_t));
			}
		}
		out_file.close();
		delete[] buf_in;
		delete[] buf_out;
	}
	for(int i = 0;i < N_THREADS;i++){
		rename(path(filename(i)+".tmp"), path(filename(i)));
	}
}

template<class T>
void coin::GloveVLBL<T>::load_glove_vector(const std::string& glove_file, const bool skip_id){
	ifstream glove_in(glove_file, ios::in);
	if(!glove_in){
		cerr << glove_file << " not found!!" << endl;
	}
	string line;
	int windex = 0;
	if(skip_id){
		getline(glove_in, line); // skip header
	}
	while(glove_in && getline(glove_in, line)){
		int dim = 0;
		istringstream iss(line);
		string word;
		iss >> word;
		vector<double> embed;
		double bias = 0.;
		if(skip_id){
			int id; double gram;
			iss >> id >> bias >> gram;
		}
		while(iss){
			double v;
			iss >> v;
			if(iss){
				embed.push_back(v);
				dim++;
			}
		}
		dim_ = dim;
		preprocess(word);
		if(word_indices_.find(word) == word_indices_.end()){
			word_indices_[word] = windex;
			word_strings_.push_back(word);
			Vector embeddings(dim);
			for(int i = 0;i < dim;i++){
				embeddings[i] = embed[i];
			}
			words_.push_back(new T(windex, embeddings, bias));
			windex++;
		}else{
			//TODO: lowercase
		}
	}
	cerr << windex << " words are loaded." << endl;
}

template<class T>
void coin::GloveVLBL<T>::save_model(const string& model_file){
	ofstream os(model_file, ios::binary | ios::out | ios::trunc);
	os.write((char *)&niters_, sizeof(int));
	os.write((char *)&dim_, sizeof(int));
	int size = words_.size();
	os.write((char *)&size, sizeof(int));
	for(T* word:words_){
		word->save(os);
	}
	assert(words_.size() == word_indices_.size());
	os.close();
}

template<class T>
void coin::GloveVLBL<T>::save_text_model(const string& model_file){
	ofstream os(model_file, ios::out | ios::trunc);
	os << dim_ << " " << words_.size() << endl;
	for(T* word:words_){
        os << word_strings_[word->id()];
		word->save_text(os);
        os << endl;
	}
	assert(words_.size() == word_indices_.size());
	os.close();
}


template<class T>
void coin::GloveVLBL<T>::load_model(const string& model_file){
	ifstream is(model_file, ios::binary | ios::in );
	if(!is){
		cerr << "error in opening model file: " << model_file << endl; 
		exit(0);
	}
	is.read((char *)&niters_, sizeof(int));
	is.read((char *)&dim_, sizeof(int));
	for(T* word:words_){
		delete word;
	}
	int size = 0;
	is.read((char *)&size, sizeof(int));
	words_.resize(size);
	for(int i = 0;i < size;i++){
		T* word = new T(dim_, is);
		words_[word->id()] = word;
	}
	is.close();
}

template<class T>
coin::GloveVLBL<T>::~GloveVLBL(){
	for(T* word:words_){
		delete word;
	}
	if(is_directory(tmpdir_)){
		remove_all(tmpdir_);
	}
}

template<class T>
void coin::GloveVLBL<T>::similar_words(const string& word, int limit, bool reverse){
	if(word_indices_.find(word) == word_indices_.end()){
		cout << "no word vector for " << word << endl;
    int i = 0;
    while(true){
      if(++i > limit)break;
      cout << "unknown" << "\t" << 0.0000 << endl;
    }
		return;
	}
	T *w = words_[word_indices_[word]];
	list<pair<int, double>> sim_words;
	for(T* comp:words_){
		if(w == comp)continue;
    if(reverse){
      sim_words.push_back(pair<int, double>(comp->id(), -w->similarity(*comp)));
    }else{
      sim_words.push_back(pair<int, double>(comp->id(), w->similarity(*comp)));
    }
	}
	sim_words.sort(more<int>);
	int i = 0;
	for(pair<int, double> sim_word:sim_words){
		if(++i > limit)break;
    if(reverse){
      cout << word_strings_[sim_word.first] << "\t" << -sim_word.second << endl;
    }else{
      cout << word_strings_[sim_word.first] << "\t" << sim_word.second << endl;
    }
	}
}

template<class T>
void coin::GloveVLBL<T>::antonym_score(const std::string& file){
	ifstream in(file.c_str());
	if(!in){
		cerr << file << " not found!!" << endl;
		return;
	}
	load_dict();
	string line;
	string target, word, ans_word;
	int w1idx, w2idx;
	T* w1;
	T* w2;
	int correct = 0, total = 0, answer = 0;
	while(in && getline(in, line)){
		total++;
		preprocess(line);
		istringstream iss(line);
		iss >> target;
		target = target.substr(0, target.size() - 1);
		if(word_indices_.find(target) != word_indices_.end()){
			w1idx = word_indices_[target];
			w1 = words_[w1idx];
		}else{
			w1idx = -1;
			continue;
		}
		ans_word = "";
		double min_sim = 1000000.;
		bool skip = false;
		for(int i = 0;i < 5;++i){
			iss >> word;
			double sim = 100.;
			if(word_indices_.find(word) != word_indices_.end()){
				w2idx = word_indices_[word];
				w2 = words_[w2idx];
				sim = w1->similarity(*w2);
			}else{
				w2idx = -1;
				skip = true;
			}
			if(sim < min_sim){
				min_sim = sim;
				ans_word = word;
			} 
			if(includes_pair(antonyms_, w1idx, w2idx)){
				cerr << target << ":" << word << ":" << sim << "*" << endl;
			}else{
				cerr << target << ":" << word << ":" << sim << endl;
			}
		}
		iss >> word;
		assert(word == "::");
		iss >> word;
		cerr << "====" << target << ":" << word << " => " << ans_word << endl;
		if(skip){
			continue;
		}
		if(ans_word == word){
			correct++;
		}
		answer++;
	}
	double p = correct/(double)answer;
	double r = correct/(double)total;
	double f = 2. * p * r / (p + r);
	cerr << "P/R/F: " << p << " / " << r << " / " << f << endl;
}

template<class T>
void coin::GloveVLBL<T>::conv(const string& gzfile){
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
		Vector v(Vector::Zero(dim_));
		while(iss >> word){
			if(word_indices_.find(word) != word_indices_.end()){
				v += words_[word_indices_[word]]->embeddings();
			}
		}
		for(int i = 0;i < dim_;++i){
			cout << v[i] << " ";
		}
		cout << endl;
	}
	in.close();
}

template<class T>
void coin::GloveVLBL<T>::iterate(){
	double cost = 0.;
	long nupdates = 0;
	Timer t;
	omp_lock_t lock[nwords_];
	for (int i = 0;i < nwords_;++i){
	  omp_init_lock(&(lock[i]));
	}

#ifdef _OPENMP
#pragma omp parallel num_threads(N_THREADS) reduction(+:nupdates) reduction(+:cost) default(shared)
#endif
	{
		//unsigned long long next_random = niters_ * omp_get_thread_num();
		char *buf = new char[BUFSIZE];
		ifstream in_file(filename(omp_get_thread_num()), ios::in | ios::binary);
		in_file.rdbuf()->pubsetbuf(buf, BUFSIZE);		
		while(true){
			int32_t w1idx, w2idx;
			count_t count, neg_count;
			in_file.read((char *)&w1idx, sizeof(int32_t));
			if(in_file.eof())break;
			in_file.read((char *)&w2idx, sizeof(int32_t));
			in_file.read((char *)&count, sizeof(count_t));
			in_file.read((char *)&neg_count, sizeof(count_t));	

			// if(count < 1. && count > 0.){
			// 	next_random = next_random * (unsigned long long)25214903917 + 11;
			// 	count = ((next_random & 0xFFFF) / (count_t)65536 < count) ? 1. : 0.;
			// }			
			// if(neg_count < 1. && neg_count > 0.){
			// 	next_random = next_random * (unsigned long long)25214903917 + 11;
			// 	neg_count = ((next_random & 0xFFFF) / (count_t)65536 < neg_count) ? 1. : 0.;
			// }		   

			T* w1 = words_[w1idx];			  
			T* w2 = words_[w2idx];			  

			// ivlbl
			double q = sigmoid(w1->t_embeddings().dot(w2->c_embeddings()) + w2->c_bias());
			//double q = sigmoid(w1->t_embeddings().dot(w2->c_embeddings())); // no bias
			double diff = count * (1. - q) - neg_count * q;
			if(q > 0. && q < 1.){
				cost -= count * log(q) + neg_count * log(1. - q);
			}
			if(diff == 0.){
			    continue;
			}
			Vector w2update(w1->t_embeddings());

			omp_set_lock(&(lock[w1idx]));
			w1->update_t(w2->c_embeddings(), 0., diff);
			omp_unset_lock(&(lock[w1idx]));
			
			omp_set_lock(&(lock[w2idx]));
			w2->update_c(w2update, 1., diff);
			//w2->update_c(w2update, 0., diff); // no bias 
			omp_unset_lock(&(lock[w2idx]));
			// #pragma omp critical
			// {
			// 	cerr << q << ":" << count << ":" << neg_count << ":" << sigmoid(w1->t_embeddings().dot(w2->c_embeddings()) + w2->c_bias()) << endl;
			//}
			nupdates++;
		}
		in_file.close();
		delete[] buf;
	}
	niters_++;
	//reset_history();

    for(int i = 0; i < nwords_;++i){
        omp_destroy_lock(&(lock[i]));
	}
	double time = (double)t.seconds();
	cerr << "cost: " << cost << ", updates: " << nupdates << ", time: " << time << " seconds, speed:" << nupdates / (N_THREADS * time) << " [words/threads/sec]" << endl;
}


template<class T>
void coin::GloveVLBL<T>::reset_history(){
	for(T* w:words_){
		w->reset_history();
	}
}


template<class T>
T coin::GloveVLBL<T>::calc_mean_word(){
	T mean(dim_);
	for(T* word:words_){
		mean.add(*word);
	}
	mean.mult(1. / words_.size());
	return mean;
}

coin::SymmetricWord::SymmetricWord(const SymmetricWord& word):
	bias_(word.bias()), grad_(word.grad_), embeddings_(word.embeddings()){};

coin::SymmetricWord& coin::SymmetricWord::operator=(const SymmetricWord& word) {
	bias_ = word.bias();
	grad_ = word.grad_;
	embeddings_ = word.embeddings();
	return *this;
}

coin::SymmetricWord::SymmetricWord(int dim):
	Word(0), bias_(0.), grad_(0.), embeddings_(Vector::Zero(dim)){
}

coin::SymmetricWord::SymmetricWord(int id, int dim, function<double()> &rand):
	Word(id), bias_(0.), grad_(0.), embeddings_(Vector(dim)){
	for(int i = 0;i < dim;++i){
		embeddings_[i] = rand();
	}
	assert(dim == embeddings_.size());
}

coin::SymmetricWord::SymmetricWord(int id, const Vector &v, double bias):
	Word(id), bias_(bias), grad_(0.), embeddings_(v){
}


inline void coin::SymmetricWord::update(const Vector& vector_grad, const double bias_grad, const double scale){
	grad_ += (vector_grad.squaredNorm() + bias_grad * bias_grad) / (vector_grad.size() + 1.) * scale * scale;
 	double alpha = RHO * scale / sqrt(grad_);
 	embeddings_ += alpha * vector_grad;
 	bias_ += alpha * bias_grad;
}

void coin::SymmetricWord::save(ostream& os){
	os.write((char *)&id_, sizeof(int));
	os.write((char *)&bias_, sizeof(double));
	os.write((char *)&grad_, sizeof(double));
	for(int i = 0;i < embeddings_.size();++i){
		os.write((char*)&embeddings_[i], sizeof(double));
	}
}

void coin::SymmetricWord::save_text(ostream& os){
	for(int i = 0;i < embeddings_.size();++i){
		os << " " << embeddings_[i];
	}
    os << " " << bias_;
}

void coin::SymmetricWord::load(istream& is){
	is.read((char *)&id_, sizeof(int));
	is.read((char *)&bias_, sizeof(double));
	is.read((char *)&grad_, sizeof(double));
	for(int i = 0;i < embeddings_.size();++i){
		is.read((char*)&embeddings_[i], sizeof(double));
	}
}

void coin::SymmetricWord::add(const SymmetricWord& word){
	embeddings_ += word.embeddings();
	bias_ += word.bias();
}

void coin::SymmetricWord::mult(const double mult){
	embeddings_ *= mult;
	bias_ *= mult;
}

// explicit instantiation of a template
namespace coin{
	template class GloveVLBL<SymmetricWord>;
};

