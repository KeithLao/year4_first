#ifndef GRADIENT_COMPUTE_H
#define GRADIENT_COMPUTE_H

#include <vector>
#include "base_matrix.h"
using namespace std;

void Gradient_compute(string algorithm, const vector<double>& received, const vector<int>& codeword, 
                      const vector<int>& syndrome, const BaseMatrix& H, vector<double>& gradient, 
                      const vector<double>& noise, const vector<double>& momentum, int iter);

#endif