#include "performance.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

void Performance::addResult(string algorithm, double snr, double fer, double ber) {
    // Check if algorithm already exists
    int algorithm_index = -1;
    for (int i = 0; i < algorithms.size(); i++) {
        if (algorithms[i] == algorithm) {
            algorithm_index = i;
            break;
        }
    }
    
    // If not exists, create new algorithm entry
    if (algorithm_index == -1) {
        algorithms.push_back(algorithm);
        algorithm_data.push_back(vector<PerformanceData>());
        algorithm_index = algorithms.size() - 1;
    }
    
    // Add data
    PerformanceData new_data;
    new_data.snr = snr;
    new_data.fer = fer;
    new_data.ber = ber;
    algorithm_data[algorithm_index].push_back(new_data);
}

void Performance::saveToFile(string code_type) {
    for (int i = 0; i < algorithms.size(); i++) {
        // Remove .txt extension if present
        string clean_code_type = code_type;
        if (clean_code_type.length() >= 4 && 
            clean_code_type.substr(clean_code_type.length() - 4) == ".txt") {
            clean_code_type = clean_code_type.substr(0, clean_code_type.length() - 4);
        }
        
        string filename = "csv/" + algorithms[i] + "_" + clean_code_type + ".csv";
        ofstream file(filename);
        
        // Write CSV header
        file << "SNR,FER,BER" << endl;
        
        // Write data
        for (int j = 0; j < algorithm_data[i].size(); j++) {
            file << fixed << setprecision(2) << algorithm_data[i][j].snr << ","
                 << scientific << setprecision(6) << algorithm_data[i][j].fer << ","
                 << scientific << setprecision(6) << algorithm_data[i][j].ber << endl;
        }
        
        file.close();
        cout << "Algorithm " << algorithms[i] << " results saved to: " << filename << endl;
    }
}