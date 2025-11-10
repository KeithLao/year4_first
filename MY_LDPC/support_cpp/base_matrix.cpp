#include "base_matrix.h"
#include <fstream>
#include <iostream>
using namespace std;

BaseMatrix::BaseMatrix(const string& filename) {
    ifstream file(filename);
    
    // Read basic information
    file >> n >> m;
    
    // Read degree information
    int degree_max[2];
    file >> degree_max[0] >> degree_max[1];
    
    degree_c.resize(n);
    degree_r.resize(m);
    for (int i = 0; i < n; i++) {
        file >> degree_c[i];
    }
    for (int i = 0; i < m; i++) {
        file >> degree_r[i];
    }
    
    // Initialize Tanner graph
    tanner_graph_j.resize(n);
    tanner_graph_i.resize(m);
    
    // Read connection relationships
    for (int i = 0; i < n; i++) {
        tanner_graph_j[i].resize(degree_c[i]);
        for (int j = 0; j < degree_c[i]; j++) {
            file >> tanner_graph_j[i][j];
            // Convert from 1-based to 0-based indexing
            tanner_graph_j[i][j]--;
        }
    }
    for (int i = 0; i < m; i++) {
        tanner_graph_i[i].resize(degree_r[i]);
        for (int j = 0; j < degree_r[i]; j++) {
            file >> tanner_graph_i[i][j];
            // Convert from 1-based to 0-based indexing
            tanner_graph_i[i][j]--;
        }
    }
    
    file.close();
    cout << "Loaded LDPC matrix: " << n << " variable nodes, " << m << " check nodes" << endl;
}