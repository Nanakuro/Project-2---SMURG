#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:43:19 2019

@author: Minh Nguyen
"""

from matplotlib import pyplot as plt
from matplotlib import mlab
import numpy as np
import itertools
import stats
'''
grid_config = "grid_config.txt"
grid_energy = "grid_energy.txt"

configs_bin, configs = [], []
bits = 0

probs = []
with open(grid_config) as f, open(grid_energy) as g:
    for c in f:
        c_list = c.strip().split()
        bits = len(c_list[0])
        c_num = int(c_list[0],2)
        configs_bin.append(c_num)
    for p in g:
        p_list = p.strip().split()
        configs.append(int(p_list[0]))
        probs.append(float(p_list[1])*len(configs_bin))
'''


        
'''
prob_pair.sort(key=lambda x: x[0])
#prob_pair = list(prob_pair for prob_pair,_ in itertools.groupby(prob_pair))
seen = set()
prob_pair = [x for x in prob_pair if x[0] not in seen and not seen.add(x[0])]
print(prob_pair)
'''
    
#print('Total probability = ', sum(prob))
#print(sum(probs))
#print(configs_bin.count(0))
#print(configs)
    
def PlotConfigHist(configs_bin, configs, probs):
    fig = plt.figure()
    plt.title('Config Histogram')
    plt.xlabel('Configurations')
    plt.ylabel('Sample probability')
    n,bins = np.histogram(configs_bin, bins=configs)#, density=True, label='data')
    mid = 0.5*(bins[1:] + bins[:-1]) - 0.5
    menStd = np.sqrt(n)
    width = 1
    plt.bar(mid, n, width=width, align='center', color='r', yerr=menStd, label='data')
    plt.plot(configs, probs, lw=0.75, label='theoretical')
    plt.legend()
    plt.show()
    #plt.savefig('3x3_grid_histogram_errorbars', dpi=300)
    
def PlotEHist(E, E_bin, b):
    fig = plt.figure()
    plt.title(r'Energy Histogram - $\beta=%.1f$' % b)
    plt.xlabel(r'Energies')
    plt.ylabel(r'Counts')
    n,bins,_ = plt.hist(E, bins=E_bin, density=False, label='data')
    mid = 0.5*(bins[1:] + bins[:-1])
    err = np.sqrt(n)
    #width = 1
    #plt.bar(mid, n, width=width, align='center', color='r', yerr=menStd, label='data')
    E_ave, E_var, E_err, E_tau = stats.Stats(np.array(E))
    #print(E_ave, E_var, E_err, E_tau)
    #plt.errorbar(mid, n, yerr=err, fmt='none', label='error')
    plt.axvline(E_ave, color='y', label='mean')
    plt.legend()
    plt.show()
    if b > 100:
        plt.savefig('energy_hist_beta_inf.png', dpi=300)
    else:
        print()
        #plt.savefig('energy_hist_beta_%.1f.png' % b, dpi=300)
    
def PlotM2Hist(M2, M2_bin, b):
    fig = plt.figure()
    plt.title(r'$M^2$ Histogram - $\beta=%.1f$' % b)
    plt.xlabel(r'$M^2$')
    plt.ylabel(r'Counts')
    #M2_bin_norm = np.linspace(min(M2_bin), max(M2_bin), 2*len(M2_bin), endpoint=True)
    n,bins,_ = plt.hist(M2, bins=M2_bin, density=False, label='data')
    mid = 0.5*(bins[1:] + bins[:-1])
    err = np.sqrt(n)#np.array([ np.sqrt(n[i]/E.count(bins[i])) for i in range(len(n)) ])
    #width = 1
    #plt.bar(mid, n, width=width, align='center', color='r', yerr=menStd, label='data')
    M2_mean, M2_var, M2_err, M2_tau = stats.Stats(np.array(M2))
    #print(M2_mean, M2_var, M2_err, M2_tau)
    #plt.errorbar(mid, n, yerr=err, fmt='none', label='error')
    plt.axvline(M2_mean, color='y', label='mean')
    plt.legend()
    plt.show()
    if b > 100:
        plt.savefig('M2_hist_beta_inf.png', dpi=300)
    else:
        print()
        #plt.savefig('M2_hist_beta_%.1f.png' % b, dpi=300)
    

c_list, E_list, m2_list = [], [], []

beta_list = np.linspace(0,1,11, endpoint=True)
beta_list = np.append(beta_list, 1E6)
print(beta_list)
for beta in beta_list:
    grid_file = 'grid_beta_%.1f.txt' % beta
    with open(grid_file,'r') as f:
        for line in f:
            line = line.strip().split()
            E_list.append(int(line[1]))
            m2_list.append(float(line[2]))
    
    E_bins = list(set(E_list))
    E_bins.sort()
    m2_bins = list(set(m2_list[:]))
    m2_bins.sort()
    
    #E_ave = sum(E_list)/len(E_list)
    #m2_ave = sum(m2_list)/len(m2_list)
    PlotEHist(E_list, E_bins, beta)
    PlotM2Hist(m2_list, m2_bins, beta)
    break
        
        
            