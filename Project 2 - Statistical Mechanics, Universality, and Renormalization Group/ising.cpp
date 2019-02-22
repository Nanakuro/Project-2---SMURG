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
#include <sys/stat.h>
#include <unistd.h>

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
    
    string getBinary() {
        string bin = "";
        for (const auto &s : all_spins) {
            bin += s==1 ? "1" : "0";
        }
        return bin;
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

void WriteGrid(int size,
               mt19937 &m,
               string spinFile, string bondFile,
               bool file_exists=false) {

    if (file_exists) {
        ofstream sp_file, b_file;
        sp_file.open(spinFile, fstream::trunc);
        b_file.open(bondFile, fstream::trunc);

        uniform_int_distribution<int> int_dist(0,1);
        uniform_real_distribution<double> r_dist(-2.0,2.0);
    
        for (int i=0; i < size*size; i++) {
            int node = i+1;
            double h = 0.0, J = 1.0;
            int spin = 1;
            //int rand = int_dist(m);
            //int spin = rand==0 ? -1 : 1;
            //h = r_dist(m);
            sp_file << node << " " << spin << " " << h << endl;
            //J = r_dist(m);
            int left = node - 1,
                right = node + 1,
                up = node - size,
                down = node + size;
            
            left % size == 0 && (left += size);
            right % size == 1 && (right -= size);
            up < 1 && (up += size*size);
            down > size*size && (down -= size*size);

            if (left+1-size == node) {
                b_file << node << " " << left << " " << J << endl;// J = r_dist(m);
            }
            if (right > node) {
                b_file << node << " " << right << " " << J << endl;// J = r_dist(m);
            }
            if (down > node) {
                b_file << node << " " << down << " " << J << endl;// J = r_dist(m);
            }
            if (up+size-size*size == node) {
                b_file << node << " " << up << " " << J << endl;// J = r_dist(m);
            }
        }
        sp_file.close();
        b_file.close();
    }
}

inline bool fileExists (const string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

double calc_alpha(double delta_E, double beta) {
    double alph = exp(-delta_E*beta);
    return alph;
}

int main(int argc, char** argv) {
    // Random device
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> real_dist(0,1);
    
    string gridSpin="grid_spin.txt", gridBond="grid_bond.txt";
    //bool not_file = (!fileExists(gridSpin) && !fileExists(gridBond));
    WriteGrid(3, mt, gridSpin, gridBond,true);//, not_file);
    
    //string spin_file="spin.txt", bond_file="bonds.txt";
    string spin_file=gridSpin, bond_file=gridBond;
    
    Ising myState (spin_file,bond_file);
    myState.printNeighbors();
    cout << "E = " << myState.getE() << endl;
    uniform_int_distribution<int> random_node(1,myState.getNumSpins());
    int num_sweeps = 10000, start_time = num_sweeps/4;
    int sweep = 10000;
    double beta = 1/2.0;


    ofstream E_file;
    string E_file_name="grid_config.txt";
    E_file.open(E_file_name, fstream::trunc);
    for (int i=0; i<num_sweeps; i++) {
        for (int j=0; j<sweep; j++) {
            int node_flip = random_node(mt);
            double alpha = calc_alpha(myState.deltaE(node_flip), beta);
            double rand_var = real_dist(mt);
            if (rand_var < alpha) {
                myState.flipSpin(node_flip);
            }
        }
        if (i+1 >= start_time) {//} && (i+1) % 10 == 0) {
            E_file << myState.getBinary() << endl;
        }
    }
    E_file.close();
    
    return 0;
}
