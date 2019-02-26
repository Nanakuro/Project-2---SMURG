#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:43:19 2019

@author: Minh Nguyen
"""

from matplotlib import pyplot as plt
import numpy as np
import itertools

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
prob_pair.sort(key=lambda x: x[0])
#prob_pair = list(prob_pair for prob_pair,_ in itertools.groupby(prob_pair))
seen = set()
prob_pair = [x for x in prob_pair if x[0] not in seen and not seen.add(x[0])]
print(prob_pair)
'''
    
#print('Total probability = ', sum(prob))
print(sum(probs))
print(configs_bin.count(0))
#print(configs)

num_bins = 2**bits

fig = plt.figure()
plt.title('Config Histogram')
plt.xlabel('Configurations')
plt.ylabel('Sample probability')
n,bins = np.histogram(configs_bin, bins=configs)#, density=True, label='data')
mid = 0.5*(bins[1:] + bins[:-1]) - 0.5
menStd = np.sqrt(n)
width = 1
print(type(n))
print(type(bins))
#plt.bar(mid, n, width=width, align='center', color='r', yerr=menStd, label='data')
plt.plot(configs, probs, lw=0.75, label='theoretical')
plt.legend()
plt.show()
#plt.savefig('3x3_grid_histogram_errorbars', dpi=300)