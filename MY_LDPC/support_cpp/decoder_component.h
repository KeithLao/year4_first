#ifndef DECODER_COMPONENT_H
#define DECODER_COMPONENT_H

#include <vector>
#include "base_matrix.h"
using namespace std;

bool syndrome_updateandcheck(const vector<int>& codeword, vector<int>& syndrome, const BaseMatrix& H);
void bit_flipping(const vector<int>& flip_labeled, vector<int>& codeword);
double cal_total_grad(const vector<double>& received, const vector<int>& codeword, const vector<int>& syndrome);
void Multi_labeling(const vector<double>& gradient, vector<int>& flip_labeled);
void Single_labeling(const vector<double>& gradient, vector<int>& flip_labeled);
void Noise_change(std::vector<double>& noise, int seed);
void Momentum_change(const vector<int>& flip_labeled, int size, vector<double>& momentum);
#endif