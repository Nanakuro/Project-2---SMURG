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
#include <cmath>

using namespace std;

class Ising {
    vector<vector<pair<int,double>>> neighbors;
public:
    vector<int> all_spins;
    vector<double> magfield;
    Ising(string, string);
    vector<vector<pair<int,double>>> getNeighbors() { return neighbors; };
    
    int getNumSpins() { return (int)all_spins.size(); }
    void flipSpin(int spin_flip) { all_spins[spin_flip-1] *= -1; }
    
    double getE() {
        double E=0.0;
        for (int i=0; i < neighbors.size(); i++) {
            int node = i+1, spin = all_spins[i];
            double h = magfield[i];
            for (int j=0; j < neighbors[i].size(); j++) {
                int neighbor = neighbors[i][j].first;
                if (node < neighbor) {
                    double J = neighbors[i][j].second;
                    int spin_neighbor = all_spins[neighbor-1];
                    //cout << "Node " << node << "\t neighbor=" << neighbor << "\t J x spin x spin_neighbor = " << J*spin * spin_neighbor << endl;
                    E += J*spin * spin_neighbor;
                }
            }
            //cout << "Node " << node << "\t spin x h = " << h*spin << endl;
            E += h*spin;
        }
        return E;
    }
    
    double deltaE(int spin_flip) {
        double deltaE = 0.0;
        
        // Complicated code
        flipSpin(spin_flip);
        int idx = spin_flip-1;
        deltaE += 2*all_spins[idx] * magfield[idx];

        for (int i=0; i < neighbors[idx].size(); i++) {
            int neighbor = neighbors[idx][i].first;
            if (spin_flip < neighbor) {
                double J = neighbors[idx][i].second;
                int spin_neighbor = all_spins[neighbor-1];
                deltaE += 2*J * all_spins[idx] * spin_neighbor;
            }
        }
        flipSpin(spin_flip);
        
        // Simple code
//        double E_beg = getE();
//        flipSpin(spin_flip);
//        double E_fin = getE();
//        flipSpin(spin_flip);
//        deltaE = E_fin - E_beg;

        return deltaE;
    }
    
    
    void printNeighbors(){
        for (int i=0; i < neighbors.size(); i++) {
            cout << "Node " << i+1 << ":";
            for (int j=0; j < neighbors[i].size(); j++) {
                pair<int,double> neighbor_pair = neighbors[i][j];
                cout << " (" << neighbor_pair.first << ", " << neighbor_pair.second << ")";
            }
            cout << endl;
        }
    }
};

Ising::Ising(string spinFile, string bondsFile) {
    ifstream spin_inp(spinFile);
    ifstream bonds_inp(bondsFile);
    int node, sp, neighbor;
    int line = 0;
    double h, J;
    while (spin_inp >> node >> sp >> h) {
        all_spins.push_back(sp);
        magfield.push_back(h);
        line++;
    }
    neighbors.resize(line);
    while (bonds_inp >> node >> neighbor >> J) {
        if (node >= 1) {
            pair<int,double> neighbor_pair = make_pair(neighbor, J);
            neighbors[node-1].push_back(neighbor_pair);
        }
        if (neighbor >= 1) {
            pair<int,double> neighbor_pair = make_pair(node, J);
            neighbors[neighbor-1].push_back(neighbor_pair);
        }
    }
}

void WriteGrid(int size, mt19937 &m, string spinFile, string bondFile) {
    ofstream sp_file, b_file;
    sp_file.open(spinFile, fstream::trunc);
    b_file.open(bondFile, fstream::trunc);

    uniform_int_distribution<int> int_dist(0,1);
    uniform_real_distribution<double> r_dist(-2.0,2.0);
    
    for (int i=0; i < size*size; i++) {
        int node = i+1;
        double h = 1.0, J = 1.0;
        int rand = int_dist(m);
        int spin = rand==0 ? -1 : 1;
        //h = r_dist(m);
        sp_file << node << " " << spin << " " << h << endl;
        //J = r_dist(m);
        if (node < size*size) {
            if (node % size != 0) {
                b_file << node << " " << node+1 << " " << J << endl;
                //J = r_dist(m);
            }
            if (node < size*size+1 - size) {
                b_file << node << " " << node+size << " " << J << endl;
                //J = r_dist(m);
            }
        }
    }
    sp_file.close();
    b_file.close();
}

int main(int argc, char** argv) {
    // Random device
    random_device rd;
    mt19937 mt(rd());
    
    string gridSpin="grid_spin.txt", gridBond="grid_bond.txt";
    WriteGrid(3,mt,gridSpin,gridBond);
    
    //string spin_file="spin.txt", bond_file="bonds.txt";
    string spin_file=gridSpin, bond_file=gridBond;
    
    Ising myState (spin_file,bond_file);
    myState.printNeighbors();
    cout << "E = " << myState.getE() << endl;
    
    return 0;
}
