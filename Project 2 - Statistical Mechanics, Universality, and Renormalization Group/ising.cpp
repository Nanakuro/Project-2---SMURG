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
#include <sstream>
#include <iomanip>

#include "input.h"

using namespace std;

const int num_sweeps = 10000;
const int start_time = 10;

string intToBin(unsigned int n) {
    string str = "";
    if (n / 2 != 0) {
        str += intToBin(n / 2);
    }
    str += to_string(n % 2);
    return str;
}

inline bool fileExists (const string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

string zfill(unsigned int z, string &s) {
    s = string(z - s.length(),'0') + s;
    return s;
}

string dtos(double db, int prec) {
    stringstream stream;
    stream << fixed << setprecision(prec) << db;
    string s = stream.str();
    return s;
}



class Ising {
    vector<vector<pair<int,double>>> neighbors;
    vector<int> original_spins;
    double Z;
public:
    vector<int> all_spins;
    vector<double> magfield;
    double beta;
    Ising(string, string, double);
    vector<vector<pair<int,double>>> getNeighbors() { return neighbors; };
    
    int getNumSpins() {
        return (int)all_spins.size();
    }
    
    void flipSpin(int spin_flip) {
        all_spins[spin_flip-1] *= -1;
    }
    
    double getE() {
        double E=0.0;
        for (int i=0; i < neighbors.size(); ++i) {
            int node = i+1, spin = all_spins[i];
            double h = magfield[i];
            for (int j=0; j < neighbors[i].size(); ++j) {
                int neighbor = neighbors[i][j].first;
                if (node < neighbor) {
                    double J = neighbors[i][j].second;
                    int spin_neighbor = all_spins[neighbor-1];
                    //cout << "Node " << node << "\t neighbor=" << neighbor << "\t J x spin x spin_neighbor = " << J*spin * spin_neighbor << endl;
                    E += -J * spin * spin_neighbor;
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
        //flipSpin(spin_flip);
        int idx = spin_flip-1;
        deltaE += -2*all_spins[idx] * magfield[idx];

        for (int i=0; i < neighbors[idx].size(); ++i) {
            int neighbor = neighbors[idx][i].first;
            double J = neighbors[idx][i].second;
            int spin_neighbor = all_spins[neighbor-1];
            deltaE += 2*J * all_spins[idx] * spin_neighbor;
        }
        
        // Simple code
//        double E_beg = getE();
//        flipSpin(spin_flip);
//        double E_fin = getE();
//        deltaE = E_fin - E_beg;
        
        //flipSpin(spin_flip);
        return deltaE;
    }
    
    void setZ(double b) {
        beta = b;
        const vector<int> temp_spins = all_spins;
        const int bit_length = (int)all_spins.size();
        string binary = string(bit_length,'0');
        Z = 0.0;
        int decimal = stoi(binary, nullptr, 2);
        while (decimal < (int) pow(2.0,bit_length)) {
            binary = intToBin(decimal);
            binary = string(bit_length - binary.length(),'0') + binary;
            binToSpins(binary);
            Z += exp(-beta*getE());
            ++decimal;
        };
        all_spins = temp_spins;
    }
    
    double getZ() {
        return Z;
    }
    
    double getAlpha(int spin_flip) {
        return exp(-beta*deltaE(spin_flip));
    }
    
    double getProb() {
        return exp(-beta*getE())/Z;
    }
    
    void binToSpins(string binary) {
        for (string::size_type i=0; i<binary.size(); ++i) {
            all_spins[i] = binary[i]=='1' ? 1 : -1;
        }
    }
    
    string getBinary() {
        string bin = "";
        for (const auto &s : all_spins) {
            bin += s==1 ? "1" : "0";
        }
        return bin;
    }
    
    void reset() {
        all_spins = original_spins;
    }
    
    void printNeighbors(){
        for (int i=0; i < neighbors.size(); ++i) {
            cout << "Node " << i+1 << ":";
            for (int j=0; j < neighbors[i].size(); ++j) {
                pair<int,double> neighbor_pair = neighbors[i][j];
                cout << " (" << neighbor_pair.first << ", " << neighbor_pair.second << ")";
            }
            cout << endl;
        }
    }
};

Ising::Ising(string spinFile, string bondsFile, double b) {
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
    original_spins = all_spins;
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
    setZ(b);
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
    
        for (int i=0; i < size*size; ++i) {
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

void writePc(string E_file_name, Ising &state) {
    ofstream E_file(E_file_name, fstream::trunc);
    vector<int> temp_vec = state.all_spins;
    string bin_string;
    for (int i=0; i<(double) pow(2.0,state.getNumSpins()); ++i) {
        bin_string = intToBin(i);
        state.binToSpins(zfill((int)state.getNumSpins(),bin_string));
        E_file << i << " " << state.getProb() << endl;
    }
    cout << state.getNumSpins() << endl;
    double config_num = pow(2.0, state.getNumSpins());
    cout << config_num << endl;
    state.all_spins = temp_vec;
    E_file.close();
}

void writeMCMC(mt19937 &m, Ising &state, int sweep, string out_file_name) {
    uniform_real_distribution<double> real_dist(0,1);
    uniform_int_distribution<int> random_node(1,state.getNumSpins());
    string beta_string = dtos(state.beta, 1);
    string out_file = out_file_name + "_" + beta_string + ".txt";
    ofstream outFile(out_file_name, fstream::trunc);
    for (int n_swp=0; n_swp<num_sweeps; ++n_swp) {
        for (int swp=0; swp<sweep; ++swp) {
            int node_flip = random_node(m);
            double alpha = state.getAlpha(node_flip);
            double rand_var = real_dist(m);
            if (rand_var < alpha) {
                state.flipSpin(node_flip);
            }
        }
        if (n_swp+1 >= start_time) {
            outFile << state.getBinary() << endl;
        }
    }
    outFile.close();
}

int main(int argc, char** argv) {
    // Random device
    random_device rd;
    mt19937 mt(rd());
    
    // Getting variables from myInputFile
    InputClass input;
    string myInputFileName = "myInputFile.txt";
    ifstream myInputFile(myInputFileName);
    input.Read(myInputFile);
    
    //double beta = input.toDouble(input.GetVariable("beta"));
    int Lx = input.toInteger(input.GetVariable("Lx"));
    int grid_size = Lx,
        N = grid_size*grid_size;
    string out_file_name = input.GetVariable("outFile"),
           spin_file_name = input.GetVariable("spinsFile"),
           bond_file_name = input.GetVariable("bondsFile");
    
    // Arbitrary variables
    int sweep = N;
    vector<double> beta_list {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0E6};

    for (const auto &beta : beta_list) {
        // Write initial configuration to spins and bonds files.
        //Currently set to s_i=1 and J_ij=1 for all i,j.
        //bool not_file = (!fileExists(spin_file_name) && !fileExists(bond_file_name));
        WriteGrid(grid_size, mt, spin_file_name, bond_file_name,true);//, not_file);

        // Making state object from files and beta.
        Ising myState (spin_file_name, bond_file_name, beta);
        //myState.printNeighbors();
        cout << "E = " << myState.getE() << endl;

        // MARKOV CHAIN MONTE CARLO:
        writeMCMC(mt, myState, sweep, out_file_name);

        // Writing p(c) to file
        //writePc(E_file_name, myState);
    }
    
    return 0;
}
