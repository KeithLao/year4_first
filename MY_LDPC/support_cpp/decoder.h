#ifndef DECODER_H
#define DECODER_H

#include <vector>
#include "base_matrix.h"
#include "performance.h"
using namespace std;

class Decoder {
public:
    BaseMatrix matrix;
    
    Decoder(BaseMatrix& H);
    ~Decoder() = default;

    vector<int> GDBF(string algorithm, const vector<double>& received, const vector<int>& hard_received, int max_iter, Performance& perf, vector<double>& noise); 
    // Input: channel initial value received, hard decision hard_received, and max iterations; 
    // Output: decoding success or result after max_iter
};

#endif