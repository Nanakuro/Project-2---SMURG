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

using namespace std;

#include "input.h"
#include "Ising.h"

void writeTheoreticalProbC(string E_file_name, Ising &state) {
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

void writeMCMC(mt19937 &m, Ising &state, int sweep, string out_file_prefix) {
    uniform_real_distribution<double> real_dist(0,1);
    uniform_int_distribution<int> random_node(1,state.getNumSpins());
    string beta_string = dtos(state.getBeta(), 1);
    string out_file_name = out_file_prefix + "_" + beta_string + ".txt";
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
            outFile << state.getBinary() << " " << state.getE() << " " << state.getM2() << endl;
        }
    }
    outFile.close();
}

void writeTestZero(mt19937 &m, Ising &state, int sweep, string out_file_prefix) {
    string beta_string = dtos(state.getBeta(), 1);
    string out_file_name = out_file_prefix + "_" + beta_string + "_test.txt";
    ofstream outFile(out_file_name, fstream::trunc);
    for (int n_swp=0; n_swp<num_sweeps; ++n_swp) {
        state.randomize(m);
        if (n_swp+1 >= start_time) {
            outFile << state.getBinary() << " " << state.getE() << " " << state.getM2() << endl;
        }
    }
    outFile.close();
}

void writeTestInf(mt19937 &m, Ising &state, int sweep, string out_file_prefix) {
    string beta_string = dtos(state.getBeta(), 1);
    string out_file_name = out_file_prefix + "_" + beta_string + "_test.txt";
    ofstream outFile(out_file_name, fstream::trunc);
    uniform_int_distribution<int> rand_int(0,1);
    int flag = rand_int(m);
    char flag_char = flag==1 ? 1 : 0;
    string state_string = string(state.getNumSpins(), flag_char);
    state.binToSpins(state_string);
    for (int n_swp=0; n_swp<num_sweeps; ++n_swp) {
        if (n_swp+1 >= start_time) {
            outFile << state.getBinary() << " " << state.getE() << " " << state.getM2() << endl;
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
    string out_name = input.GetVariable("outFile"),
           spin_name = input.GetVariable("spinsFile") + ".txt",
           bond_name = input.GetVariable("bondsFile") + ".txt";
    
    // Arbitrary variables
    int sweep = N;
    vector<double> beta_list {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0E6};

    // Loop through beta_list
//    for (const auto &beta : beta_list) {
//        // Write initial configuration to spins and bonds files.
//        //Currently set to s_i=1 and J_ij=1 for all i,j.
//        //bool not_file = (!fileExists(spin_file_name) && !fileExists(bond_file_name));
//        WriteGrid(grid_size, mt, spin_name, bond_name,true);//, not_file);
//
//        // Making state object from files and beta.
//        Ising myState(spin_name, bond_name, beta);
//        //myState.printNeighbors();
//        cout << "E = " << myState.getE() << endl;
//
//        // MARKOV CHAIN MONTE CARLO:
//        //writeMCMC(mt, myState, sweep, out_name);
//
//        // Writing p(c) to file
//        //writeTheoreticalProbC(E_file_name, myState);
//    }
    
    // Write test sample for beta=0.0 and beta=inf
    Ising myState(spin_name, bond_name, 0.0);
    writeTestZero(mt, myState, sweep, out_name);
    
    myState.setBeta(1E6);
    writeTestInf(mt, myState, sweep, out_name);
    
    return 0;
}
