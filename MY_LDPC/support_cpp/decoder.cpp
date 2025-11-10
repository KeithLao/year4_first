#include "decoder.h"
#include "decoder_component.cpp"
#include "Gradient_compute.cpp"
#include <cmath>
#include <algorithm>
#include <vector>
using namespace std;

Decoder::Decoder(BaseMatrix& H) : matrix(H) {
}

vector<int> Decoder::GDBF(string algorithm, const vector<double>& received, const vector<int>& hard_received, int max_iter, Performance& perf, vector<double>& noise) {
    int n = matrix.n;
    int m = matrix.m;
    
    // Initialize
    vector<int> codeword = hard_received;  // Start from hard decision
    vector<int> syndrome(m, 0);
    vector<double> gradient(n, 0.0);
    vector<double> momentum(n, 0.0);

    double total_grad_now = 0.0;
    bool label_switch = false;

    // GDBF iteration
    for (int iter = 0; iter < max_iter; iter++) {
        // Calculate gradient for each variable node
        
        bool finish = syndrome_updateandcheck(codeword, syndrome, matrix);
        if (finish){
            break;
        }
        
        Gradient_compute(algorithm, received, codeword, syndrome, matrix, gradient, noise, momentum, iter);

        vector<int> flip_labeled(n, 0);
        
        double total_grad_next = cal_total_grad(received, codeword, syndrome); 
        if(total_grad_now < total_grad_next && label_switch == false){
            Multi_labeling(gradient, flip_labeled);
        }
        else{
            label_switch = true;
            Single_labeling(gradient,flip_labeled);
        }
        total_grad_now = total_grad_next;
        

        Noise_change(noise, iter);
        Momentum_change(flip_labeled, n, momentum);

        bit_flipping(flip_labeled, codeword);
    }

    return codeword;
}