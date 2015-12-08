/*
 * Utils.h
 *
 *  Created on: 2014/09/14
 *      Author: miwa
 */
#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <cstdint>
#include <unordered_map>
#include <chrono>
#include <vector>

#include <eigen3/Eigen/Dense>

namespace coin{

using Vector = Eigen::VectorXd;

// preprocess strings
extern void preprocess(std::string& str);
extern double load_word_counts(const std::string& file, std::unordered_map<std::string, int32_t> &word_indices, std::vector<std::string> &word_strings, std::vector<double> &word_counts);
extern double cosine(const Vector& l, const Vector& r);
extern std::pair<int32_t, int32_t> to_pair(int64_t i);
extern int64_t to_scalar(int32_t i1, int32_t i2);
extern double sigmoid(const double x);
extern void print_mem();

class Timer{
private:
	std::chrono::time_point<std::chrono::system_clock> start_;
public:
    Timer(): start_(std::chrono::system_clock::now()){}
	int seconds() const{
        return (int)std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start_).count();
	}
};

};

#endif /* UTILS_H_ */
