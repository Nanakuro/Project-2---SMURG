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
    string out_file_name = "grid_RG/" + out_file_prefix + "_" + beta_string + ".txt";
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
            outFile << state.getM2() << endl;
        }
    }
    outFile.close();
//    cout << "...MCMC done" << endl;
}

void writeMCMC(mt19937 &m, Ising &state, int sweep, string out_file_prefix, vector<double> beta_list) {
    uniform_real_distribution<double> real_dist(0,1);
    uniform_int_distribution<int> random_node(1,state.getNumSpins());
    string beta_string = dtos(state.getBeta(), 1);
    string out_file_name = "grid_RG/" + out_file_prefix + ".txt";
    ofstream outFile(out_file_name, fstream::trunc);
    for (const auto &b : beta_list) {
        state.reset();
        state.setBeta(b);
        outFile << state.getBeta() << " ";
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
                outFile << state.getM2() << " ";
            }
        }
        outFile << endl;
    }
    outFile.close();
    //    cout << "...MCMC done" << endl;
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

void writeCv(mt19937 &m, Ising &state, int sweep, string out_file_prefix, vector<double> beta_list) {
    uniform_real_distribution<double> real_dist(0,1);
    uniform_int_distribution<int> random_node(1,state.getNumSpins());
    string file_string = "Cv";
    string out_file_name = out_file_prefix + "_" + file_string + ".txt";
    ofstream outFile(out_file_name, fstream::trunc);
    for (const auto &beta : beta_list) {
        state.reset();
        state.setBeta(beta);
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
                outFile << state.getE() << " ";
            }
        }
        outFile << 1/state.getBeta() << endl;
    }
    outFile.close();
    cout << "...done" << endl;
}

void writeTestCG(mt19937 &m, Ising &state, int cg_scale, string out_file_prefix) {
    string out_file_name = "grid_RG/" + out_file_prefix + "_0.0_compare_cg_test.txt";
    ofstream outFile(out_file_name, fstream::trunc);
    //outFile << "Nat CG" << endl;
    for (int swp=0; swp<num_sweeps; ++swp) {
        state.randomize(m);
        state.setNewDefault();
        state.autoCG(cg_scale);
        outFile << state.getM2() << " ";
        state.randomize(m);
        outFile << state.getM2() << endl;
        state.reset();
    }
    outFile.close();
}

void writeCG(mt19937 &m, Ising &state, int sweep, string out_file_prefix, int cg_scale, vector<double> beta_list) {
    uniform_real_distribution<double> real_dist(0,1);
    uniform_int_distribution<int> random_node(1,state.getNumSpins());
//    string out_file_name = "grid_RG/" + out_file_prefix + "_cg.txt";
//    ofstream outFile(out_file_name, fstream::trunc);
    
    string beta_string = dtos(state.getBeta(), 1);
    string out_compare = "grid_RG/" + out_file_prefix + "_compare_cg.txt";
    ofstream outCompare(out_compare, fstream::trunc);
    for (const auto &b : beta_list) {
        state.reset();
        state.setBeta(b);
//        outFile << state.getBeta() << " ";
        for (int n_swp=0; n_swp<num_sweeps; ++n_swp) {
            for (int swp=0; swp<sweep; ++swp) {
                int node_flip = random_node(m);
                double alpha = state.getAlpha(node_flip);
                double rand_var = real_dist(m);
                if (rand_var < alpha) {
                    state.flipSpin(node_flip);
                }
            }
//            if (n_swp+1 >= start_time) {
//                vector<int> temp = state.all_spins;
//                state.autoCG(cg_scale);
//                outFile << state.getM2() << " ";
//                state.all_spins = temp;
//            }
        }
        if (beta_list.size() <= 10) {
            outCompare << state.getBeta() << " " << state.getBinary() << " ";
            state.autoCG(cg_scale);
            outCompare << state.getBinary() << " ";
            state.autoCG(cg_scale);
            outCompare << state.getBinary() << endl;
        }
//        outFile << endl;
    }
//    outFile.close();
    outCompare.close();
//    cout << "...CG done" << endl;
}

struct InputParams {
    int Lx;
    int cg_scale;
    string out_file_name;
    string spin_file_name;
    string bond_file_name;
};

InputParams getInput(string inputFileName) {
    InputParams input_params;
    
    InputClass input;
    ifstream inputFile(inputFileName);
    input.Read(inputFile);
    
    input_params.Lx = input.toInteger(input.GetVariable("Lx"));
    input_params.cg_scale = input.toInteger(input.GetVariable("cg_size"));
    input_params.out_file_name = input.GetVariable("outFile");
    input_params.spin_file_name = input.GetVariable("spinsFile") + ".txt";
    input_params.bond_file_name = input.GetVariable("bondsFile") + ".txt";
    
    return input_params;
}

int main(int argc, char** argv) {
    // Random device
    random_device rd;
    mt19937 mt(rd());
    

    
    // Getting variables from myInputFile
    string myInputFileName = "myInputFileRG.txt";
    InputParams params = getInput(myInputFileName);
    int     grid_size   = params.Lx,
            cg_size     = params.cg_scale;
    string  out_name    = params.out_file_name,
            spin_name   = params.spin_file_name,
            bond_name   = params.bond_file_name;
    
    // Arbitrary variables
    int sweep = grid_size*grid_size;//, N_beta = 200;
    
    //vector<double> beta_list {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0E6};
    
//    vector<double> T_list = makeSequenceVec(1.25, 10, N_beta, true);
//    vector<double> beta_list(T_list.size());
//    for (int i=0; i<T_list.size(); ++i) {
//        beta_list[beta_list.size()-1-i] = 1/T_list[i];
//    }

    vector<double> beta_list {0.0,0.3,0.4,0.5,0.6,1.0E6};
    //vector<double> beta_list = makeSequenceVec(0.0, 1.0, 101, true);
    
    WriteGrid(grid_size, spin_name, bond_name,true);//, not_file);

    // Making state object from files and beta.
    Ising myStateCG(spin_name, bond_name, beta_list[0]);
    
    // Write the MCMC energies and betas into a file
//    writeCv(mt, myState, sweep, out_name, beta_list);

    // Write a random grid to compare coarse graining
    writeTestCG(mt, myStateCG, cg_size, out_name);

    // Write coarse grained:
//    writeCG(mt, myStateCG, sweep, out_name, cg_size, beta_list);
    
//    WriteGrid((int)round(grid_size/cg_size), spin_name, bond_name, true);//, not_file);
//    Ising myState(spin_name, bond_name, beta_list[0]);
    
    //writeMCMC(mt, myState, sweep, out_name, beta_list);
    
////////////////////////////////////////////////////////////////////////////////////////////////
//    uniform_real_distribution<double> real_dist(0,1);
//    uniform_int_distribution<int> random_node(1,myState.getNumSpins());
//    uniform_int_distribution<int> random_node_cg(1,myStateCG.getNumSpins());
//    string out_file_name    = "grid_RG/" + out_name + ".txt";
//    string out_file_name_cg = "grid_RG/" + out_name + "_cg.txt";
//    ofstream outFile    (out_file_name, fstream::trunc);
//    ofstream outFileCG  (out_file_name_cg, fstream::trunc);
//
//    for (const auto &b : beta_list) {
//        myState.reset();
//        myState.setBeta(b);
//        myStateCG.reset();
//        myStateCG.setBeta(b);
//
//        outFile     << myState.getBeta()    << " ";
//        outFileCG   << myStateCG.getBeta()  << " ";
//        for (int n_swp=0; n_swp<num_sweeps; ++n_swp) {
//            for (int swp=0; swp<sweep; ++swp) {
//
//                int node_flip       = random_node(mt);
//                int node_flip_cg    = random_node_cg(mt);
//
//                double alpha        = myState.getAlpha(node_flip);
//                double alpha_cg     = myStateCG.getAlpha(node_flip_cg);
//
//                double rand_var     = real_dist(mt);
//                double rand_var_cg  = real_dist(mt);
//
//                if (rand_var < alpha)           { myState.flipSpin(node_flip); }
//                if (rand_var_cg < alpha_cg)     { myStateCG.flipSpin(node_flip_cg); }
//            }
//            if (n_swp+1 >= start_time) {
//                vector<int> temp = myStateCG.all_spins;
//                myStateCG.autoCG(cg_size);
//                outFileCG << myStateCG.getM2() << " ";
//                myStateCG.all_spins = temp;
//
//                outFile << myState.getM2() << " ";
//            }
//        }
//        outFile << endl;
//        outFileCG << endl;
//    }
//    outFile.close();
//    outFileCG.close();
////////////////////////////////////////////////////////////////////////////////////////////////

    // Loop through beta_list
//    for (const auto &beta : beta_list) {
//        // Write initial configuration to spins and bonds files.
//        //Currently set to s_i=1 and J_ij=1 for all i,j.
//
//        //bool not_file = (!fileExists(spin_file_name) && !fileExists(bond_file_name));
//        WriteGrid(grid_size, spin_name, bond_name, true);//, not_file);
//        Ising myState(spin_name, bond_name, beta);
//
//        WriteGrid((int)round(grid_size/cg_size), spin_name, bond_name, true);//, not_file);
//        myState = Ising(spin_name, bond_name, beta);
//
//        // MARKOV CHAIN MONTE CARLO:
//        writeMCMC(mt, myState, sweep, out_name);
//
//        // Write p(c) to file:
//        //writeTheoreticalProbC(E_file_name, myState);
//
//        cout << beta << endl << endl;
//    }
    
    // Write test sample for beta=0.0 and beta=inf
//    Ising myState(spin_name, bond_name, 0.0);
//    writeTestZero(mt, myState, sweep, out_name);
//
//    myState.setBeta(1E6);
//    writeTestInf(mt, myState, sweep, out_name);
    
    return 0;
}
