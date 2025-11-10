#ifndef INI_SEQ_H
#define INI_SEQ_H

#include "base_matrix.h"
#include <vector>
using namespace std;

class IniSeq {
public:
    int n,m,k;
    
    IniSeq(int a, int b, int c);
    
    vector<int> generateRandom(int seed); // Generate random 0,1 sequence of length bits based on seed
    vector<int> Encoder(const BaseMatrix& H, const vector<int>& ini);    // Encode the ini to message
    vector<double> addAWGN(const vector<int>& input, double snr, int seed, double bitrate); 
    vector<double> generate_noise(double snr);
    // Input is 1, -1 sequence of length bits, add noise to sequence based on snr strength, 
    // Return double sequence of same length
    vector<int> hard_decision(const vector<double>& received);
};

#endif