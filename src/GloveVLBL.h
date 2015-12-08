/*
 * GloveVLBL.h
 *
 *  Created on: 2014/09/13
 *      Author: miwa
 */
#ifndef GLOVE_VLBL_H_
#define GLOVE_VLBL_H_

#include <unordered_map>
#include <vector>
#include <string>
#include <functional>
#include <set>
#include <map>
#include "Utils.h"

namespace coin{

class Word {
protected:
	static constexpr double RHO = 3.e-2;
	static constexpr double EPSILON = 1e-12;
	int id_;
public:
	Word():id_(0){}
	Word(int id):id_(id){}
	virtual ~Word(){}
	virtual void reset_history() = 0;
	virtual double bias() const = 0;
	virtual double t_bias() const = 0;
	virtual double c_bias() const = 0;
	virtual const Vector& t_embeddings() const = 0;
	virtual const Vector& c_embeddings() const = 0;
	virtual void update_t(const Vector& vector_grad, const double bias_grad, const double scale) = 0;
	virtual void update_c(const Vector& vector_grad, const double bias_grad, const double scale) = 0;
	virtual void save(std::ostream& os) = 0;
	virtual void save_text(std::ostream& os) = 0;
	virtual void load(std::istream& is) = 0;
	int id() const{
		return id_;
	}
};

class SymmetricWord final : public Word{
private:
	double bias_ = 0.;
	double grad_ = EPSILON;
	Vector embeddings_;
	void update(const Vector& vector_grad, const double bias_grad, const double scale);
public:
	SymmetricWord(int dim);
	SymmetricWord(int id, int dim, std::function<double()> &rand);
	SymmetricWord(int id, const Vector &v, const double bias = 0.);
	SymmetricWord(int dim, std::istream& is):embeddings_(Vector(dim)){
		load(is);
	};
	SymmetricWord(const SymmetricWord& word);
	SymmetricWord& operator=(const SymmetricWord& word);
	~SymmetricWord(){}
	void reset_history(){
		grad_ = EPSILON;
	}
	double bias() const{
		return bias_;
	}
	double t_bias() const {
		return bias_;
	}
	double c_bias() const {
		return bias_;
	}
	const Vector& embeddings() const {
		return embeddings_;
	}
	const Vector& t_embeddings() const {
		return embeddings_;
	}
	const Vector& c_embeddings() const {
		return embeddings_;
	}
	void update_t(const Vector& vector_grad, const double bias_grad, const double scale){
		update(vector_grad, bias_grad, scale);
	}
	void update_c(const Vector& vector_grad, const double bias_grad, const double scale){
		update(vector_grad, bias_grad, scale);
	}
	double similarity(const SymmetricWord &word) const{
		return (sigmoid(word.t_embeddings().dot(c_embeddings()) + c_bias()) + sigmoid(t_embeddings().dot(word.c_embeddings()) + word.c_bias())) / 2.;
	}
	void save(std::ostream& os);
	void save_text(std::ostream& os);
	void load(std::istream& is);
	void add(const SymmetricWord& word);
	void mult(const double mult);
};

typedef double count_t;

class UpdateEntry{
public:
	int32_t w1;
	int32_t w2;
	count_t count;
	count_t neg_count;
 UpdateEntry(int32_t w1_, int32_t w2_, count_t count_, count_t neg_count_):
	w1(w1_), w2(w2_), count(count_), neg_count(neg_count_){}
};

template<class T>
class GloveVLBL final{
private:
	static constexpr int N_THREADS = 10;
	static constexpr double K = 5.;
	static constexpr count_t MIN_COUNT = 1.e-5;
	static constexpr int FILES = N_THREADS * 4;
	static constexpr int BUFSIZE = 1024 * 1024;
	static constexpr int SHUFFLE_SIZE = (int)(128. * 1024 * 1024 / (2. * sizeof(int) + 2 * sizeof(count_t)));
	static constexpr double SUBSAMPLING = 1.e-8;
	//static constexpr double SUBSAMPLING = 1.;
	static constexpr count_t SYN_WEIGHT = 1.e2;
	static constexpr count_t ANT_WEIGHT = (1305702. / 408302.) * SYN_WEIGHT;
	double *expTable;
	std::string tmpdir_;
	std::string coocc_file_;
	int niters_ = 0;
	int dim_ = 0.;
	int nwords_ = 0;
	double total_ = 0.;
	std::unordered_map<std::string, int> word_indices_;
	std::vector<std::string> word_strings_;
	std::vector<T*> words_;
	std::vector<double> word_counts_;
	std::unordered_map<int, std::set<int> > antonyms_;
	std::unordered_map<int, std::set<int> > synonyms_;
	void merge_files();
	void load_dict(const std::string& dict, std::unordered_map<int, std::set<int> >& map);
	void load_dict();
	void load_glove_vector(const std::string& glove_file, const bool skip_id);
	bool includes_pair(std::unordered_map<int, std::set<int> >& map, int w1, int w2);
	T calc_mean_word();
	void reset_history();
	void init_table();
	std::string filename(int idx) const;
public:
	GloveVLBL(const std::string& freq_file, const std::string& coocc_file, int dim);
	GloveVLBL(const std::string& freq_file, const std::string& model_file);
	GloveVLBL(const std::string& freq_file, const std::string& coocc_file, const std::string& model_file);
	~GloveVLBL();
	void iterate();
	void conv(const std::string& gzfile);
	void antonym_score(const std::string& file);
	void similar_words(const std::string& word, int limit = 10, bool reverse=false);
	void save_model(const std::string& model_file);
	void save_text_model(const std::string& model_file);
	void load_model(const std::string& model_file);
};

}

#endif /* GLOVE_VLBL_H_ */
