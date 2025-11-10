#include "decoder_component.h"
#include <cmath>
#include <algorithm>
#include <vector>
using namespace std;

bool syndrome_updateandcheck(const vector<int>& codeword, vector<int>& syndrome, const BaseMatrix& H){
    // Update syndrome - need to recalculate all affected syndromes
    int m = H.m;
    syndrome.assign(m, 1);  // Reset syndrome
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < H.degree_r[i]; j++) {
            int var_node = H.tanner_graph_i[i][j];
            syndrome[i] *= (double)codeword[var_node];
        }
        //syndrome[i] = (syndrome[i] + 1) >> 1;
        // adding masking
    }
    
    // Check if done
    bool finish = true;
    for (int i = 0; i < m; i++) {
        if (syndrome[i] != 1) {
            finish = false;
            break;
        }
    }
    return finish;
}

void bit_flipping(const vector<int>& flip_labeled, vector<int>& codeword){
    // Flip all bits decoder choosed
    for (int idx : flip_labeled) {
        codeword[idx] = -codeword[idx];
    }    
}


double cal_total_grad(const vector<double>& received, const vector<int>& codeword, const vector<int>& syndrome){
    int n = received.size();
    int m = syndrome.size();
    double total_grad = 0.0;
    for(int idx = 0; idx < n; idx++){
        if(idx < m){
            total_grad += syndrome[idx];
        }
        total_grad += received[idx] * codeword[idx];
    }
    return total_grad;
}

void Multi_labeling(const vector<double>& gradient, vector<int>& flip_labeled){
    int n = gradient.size();
    Single_labeling(gradient, flip_labeled);

    double min_grad = gradient[flip_labeled[0]];
    flip_labeled.pop_back();

    double flip_threshold = 0.25;
    // Flip multi. bits
    for (int j = 0; j < n; j++) {
        if (gradient[j] <= (min_grad + flip_threshold)){  
            flip_labeled.push_back(j);
        }
    }    
}

void Single_labeling(const vector<double>& gradient, vector<int>& flip_labeled){
    // Find minimum gradient value
    int n = gradient.size();
    flip_labeled.clear();

    double min_gradient = gradient[0];
    int label = 0;
    for (int j = 1; j < n; j++) {
        if (gradient[j] < min_gradient) {
              min_gradient = gradient[j];
              label = j;
        }
    }
    flip_labeled.push_back(label);   
}

void Noise_change(std::vector<double>& noise, int seed) {
    std::mt19937 rng(seed);
    std::shuffle(noise.begin(), noise.end(), rng);
}

void Momentum_change(const vector<int>& flip_labeled, int size, vector<double>& momentum){
    // momentum record 1 iter
    momentum.clear();
    momentum.resize(size, 0.0);
    double one_value = 1.0;

    for(int i : flip_labeled){
        momentum[i] = one_value;
    }
}