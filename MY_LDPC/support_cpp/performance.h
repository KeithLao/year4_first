#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#include <vector>
#include <string>
using namespace std;

struct PerformanceData{
    double snr;
    double fer;
    double ber;
};

class Performance {
public:
    vector<string> algorithms;                      // Outer category: store all algorithm names
    vector<vector<PerformanceData>> algorithm_data; // Data corresponding to each algorithm

    Performance() {}

    void addResult(string algorithm, double snr, double fer, double ber); 
    // Store data result of some algorithm in class
    void saveToFile(string code_type);    
    // For each algorithm in class, generate a .csv file named "{algorithm}_{code_type}.csv", 
    // Content is vector<data> results
};

#endif