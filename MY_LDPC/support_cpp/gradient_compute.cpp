#include "gradient_compute.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
using namespace std;

void Gradient_compute(string algorithm, const vector<double>& received, const vector<int>& codeword, const vector<int>& syndrome, const BaseMatrix& H, vector<double>& gradient, const vector<double>& noise, const vector<double>& momentum, int iter){

    //double p = 0.9;
    //static thread_local std::mt19937 rng_b(std::random_device{}()); 
    //std::bernoulli_distribution dist_b(p);

    // Calculate gradient for each variable node
    for (int j = 0; j < H.n; j++) {
        // Calculate syndrome of check nodes connected to variable node j
        int sum_syndrome = 0;
        int sum_syndrome_mask = 0;
        for (int k = 0; k < H.degree_c[j]; k++) {
            int check_node = H.tanner_graph_j[j][k];
            sum_syndrome += syndrome[check_node];

            //bool mask = dist_b(rng_b) ? 1 : 0;
        }
        
        // Calculate gradient:
        if(algorithm == "GDBF"){
            gradient[j] = received[j] * codeword[j] + sum_syndrome;
        }
        else if(algorithm == "NGDBF"){
            gradient[j] = received[j] * codeword[j] + sum_syndrome + noise[j];
        }
        else if(algorithm == "UPGDBF"){
            gradient[j] = received[j] * codeword[j] + sum_syndrome + momentum[j] + noise[j];
            if(iter < 20){
                gradient[j] -= noise[j];
            }
        }
        else if(algorithm == "MethodA"){
            
        }
    }
    
}