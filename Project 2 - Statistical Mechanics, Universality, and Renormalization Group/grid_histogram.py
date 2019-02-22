#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:43:19 2019

@author: Minh Nguyen
"""

from matplotlib import pyplot as plt
import numpy as np

grid_config = "grid_config.txt"
c_list = []
with open(grid_config) as f:
    for c in f:
        c_num = int(c,2)
        c_list.append(c_num)
        

num_bins = 2**len(bin(c_list[0])[2:])
fig = plt.figure()
plt.title('Config Histogram')
plt.xlabel('Configurations')
plt.ylabel('Counts')
plt.hist(c_list, bins=num_bins)
plt.show()