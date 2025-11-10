#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <vector>
#include <string>
using namespace std;

class BaseMatrix {
public:
    int n, m;  // Number of variable nodes, number of check nodes
    vector<int> degree_c, degree_r;  // Variable node degrees, check node degrees
    vector<vector<int>> tanner_graph_j, tanner_graph_i;  // Tanner graph connections
    
    BaseMatrix(const string& filename);
    ~BaseMatrix() = default;
};

#endif