//  Compile:    g++ -O3 -std=c++11 ising.cpp -o ising
//  File name:  ising.cpp
//  Project 2 - Statistical Mechanics, Universality, and Renormalization Group
//
//  Created by Minh on 02/20/19.
//

#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <utility>

using namespace std;

class Ising {
private:
    vector<vector<pair<int,double>>> neighbors;
public:
    vector<int> spins;
    vector<double> magfield;
    Ising(string, string);
    vector<vector<pair<int,double>>> getNeighbors() { return neighbors; };
};

Ising::Ising(string spinFile, string bondsFile) {
    ifstream spin_inp(spinFile);
    ifstream bonds_inp(bondsFile);
    int idx, sp, neighbor;
    int line = 0;
    double h, J;
    while (spin_inp >> idx >> sp >> h) {
        spins.push_back(sp);
        magfield.push_back(h);
        line++;
    }
    neighbors.resize(line);
    while (bonds_inp >> idx >> neighbor >> J) {
        if (idx >= 1) {
            neighbors[idx-1].push_back(make_pair(idx, J));
        }
    }
}

int main() {
    Ising myState ("spin.txt","bonds.txt");
    for (int i=0; i < myState.getNeighbors().size(); i++){
        cout << "Node " << i+1 << ":";
        for (const auto &j : myState.getNeighbors()[i]) {
            cout << " (" << j.first << ", " << j.second << ")";
        }
        cout << endl;
    }
    return 0;
}
