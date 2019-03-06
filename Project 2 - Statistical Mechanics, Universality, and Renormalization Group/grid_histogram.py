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
 
def Gaussian(x, mu, var):
    sigma = var**(0.5)
    y = (1/(2*np.pi*sigma**2)) * np.exp(-(x-mu)**2/(2*sigma**2))
    return y
   
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
    
def PlotEHist(E, E_bin, b, test=[], save=False):
    fig = plt.figure()
    plt.xlabel(r'Energies')
    plt.ylabel(r'$P(E)$')
    n,bins = np.histogram(E, bins=E_bin, density=True)
    #err = np.sqrt(n)
    width = 4#int(round((bins[-1] - bins[0])/len(bins)))
    mid = 0.5*(bins[1:] + bins[:-1])
    n *= width
    plt.bar(mid, n, width=width, align='center', label='data')
    #print(mid)
    E_mean, E_var, E_err, E_tau = stats.Stats(np.array(E))
    #print(E_ave, E_var, E_err, E_tau)
    plt.errorbar(E_mean, max(n)/2, color='r', capsize=5, xerr=E_err, fmt='none', label='error')
    plt.axvline(E_mean, lw=0.25, color='y', label='mean')
    
    #binWidth = (E_bin[-1] - E_bin[0])//len(E_bin)
    #x_arr = np.linspace(E_bin[0], E_bin[-1], 300, endpoint=True)
    #y_arr = Gaussian(x_arr, E_mean, E_var) * len(E)**(0.5)
    #x_arr += binWidth//2
    #plt.plot(x_arr, y_arr, color='r', label='curve')
    
    if len(test) > 0:
        PlotETest(test, E_bin, b)
    
    plt.legend()
    if b > 1E4:
        plt.title(r'Energy Histogram - $\beta=\infty$')
        save_name = 'img/energy_hist_beta_infty.png'
    else:
        #print()
        plt.title(r'Energy Histogram - $\beta=%.1f$' % b)
        save_name = 'img/energy_hist_beta_%.1f.png' % b
    
    if save:
        plt.savefig(save_name, dpi=300)
    else:    
        plt.show()
    
def PlotM2Hist(M2, M2_bin, b, test=[], save=False):
    fig = plt.figure()
    plt.xlabel(r'$M^2/N$')
    plt.ylabel(r'$P(M^2)$')
    n,bins = np.histogram(M2, bins=M2_bin, density=True)
    #err = np.sqrt(n)#np.array([ np.sqrt(n[i]/E.count(bins[i])) for i in range(len(n)) ])
    width = 0.5#int(round((bins[-1]-bins[0])/len(bins)))
    mid = 0.5*(bins[1:] + bins[:-1])
    n *= width
    plt.bar(mid, n, width=width, align='center', label='data')
    M2_mean, M2_var, M2_err, M2_tau = stats.Stats(np.array(M2))
    #print(M2_mean, M2_var, M2_err, M2_tau)
    plt.errorbar(M2_mean, max(n)/2, color='r', capsize=5, xerr=M2_err, fmt='none', label='error')
    plt.axvline(M2_mean, lw=0.25, color='y', label='mean')
    
    #binWidth = (M2_bin[-1] - M2_bin[0])//len(M2_bin)
    #x_arr = np.array(sorted(list(set(abs(np.linspace(M2_bin[0], M2_bin[-1], 300, endpoint=True))))))
    #y_arr = Gaussian(x_arr, 0, M2_var) * len(M2)**(0.5)
    #x_arr += binWidth//2
    #plt.plot(x_arr, y_arr, color='r', label='curve')
    
    if len(test) > 0:
        PlotM2Test(test, M2_bin, b)
    
    plt.legend()
    if b > 1E4:
        plt.title(r'$M^2/N$ Histogram - $\beta=\infty$')
        save_name = 'img/M2_hist_beta_infty.png'
    else:
        #print()
        plt.title(r'$M^2/N$ Histogram - $\beta=%.1f$' % b)
        save_name = 'img/M2_hist_beta_%.1f.png' % b
    if save:
        plt.savefig(save_name, dpi=300)
    else:
        plt.show()

def PlotETest(E, E_bin, b):
    n,bins = np.histogram(E, bins=E_bin, density=True)
    width = 2
    mid = 0.5*(bins[1:] + bins[:-1])
    n *= 2*width
    plt.bar(mid, n, width=width, align='center', label='test')
    E_mean, E_var, E_err, E_tau = stats.Stats(np.array(E))
    plt.axvline(E_mean, lw=0.25, label='test mean')

def PlotM2Test(M2, M2_bin, b):
    n,bins = np.histogram(M2, bins=M2_bin, density=True)
    width = 0.25
    mid = 0.5*(bins[1:] + bins[:-1])
    n *= 2*width
    plt.bar(mid, n, width=width, align='center', label='test')
    M2_mean, M2_var, M2_err, M2_tau = stats.Stats(np.array(M2))
    plt.axvline(M2_mean, lw=0.25, label='test mean')

def PlotMat(string, size, beta):
    if len(string) == size**2:
        str_split = [ c for c in string ]
        str_list = [ str_split[i:i+size] for i in range(0, len(str_split),size) ]
        str_list = np.array(str_list, int)
        fig = plt.figure()
        if beta > 1E4:
            plt.title(r'Snapshot $\beta = \infty$')
            save_name = 'img/grid_snapshot_beta_inf.png'
        else:
            plt.title(r'Snapshot $\beta = %.1f$' % beta)
            save_name = 'img/grid_snapshot_beta_%.1f.png' % beta
        plt.imshow(str_list)
        plt.savefig(save_name, dpi=300)
    else:
        print('Cannot split string into a %gx%g grid' % (size, size))
    

'''    GRAPH E & M^2 BETA = 0.0, 0.1,..., 1.0, 1000000.0    '''

beta_list = np.linspace(0,1,11, endpoint=True)
beta_list = [ '%.1f' % b for b in beta_list ]
beta_list.append('1000000.0')
#beta_list = [ str(b) for b in beta_list ]
#print(beta_list)
beta_test = ['0.0', '1000000.0']

for beta in beta_list:
    c_list, E_list, m2_list = [], [], []
    #beta = beta_list[0]
    #print(beta)
    grid_file = 'grid_beta_%s.txt' % beta
    with open(grid_file,'r') as f:
        for line in f:
            line = line.strip().split()
            c_list.append(line[0])
            E_list.append(int(line[1]))
            m2_list.append(float(line[2]))
    
    
    E_bins = np.arange(-2*len(c_list[0]), 2*len(c_list[0])+1, 4)
    m2_bins = np.arange(0, max(m2_list)+0.5, 0.5)
    #print(m2_bins)
    
    #E_bins = list(set(E_list))
    #E_bins.sort()
    #m2_bins = list(set(m2_list[:]))
    #m2_bins.sort()
    #
    #if len(E_bins) == 1:
    #    E_bins.append(E_bins[0] + 1)
    #if len(m2_bins) == 1:
    #    m2_bins.append(m2_bins[0] + 1)
    
    #E_ave = sum(E_list)/len(E_list)
    #m2_ave = sum(m2_list)/len(m2_list)
    #PlotEHist(E_list, E_bins, float(beta))
    #PlotM2Hist(m2_list, m2_bins, float(beta))
    PlotMat(c_list[-1], 27, float(beta))
    
    
    if beta in beta_test:
        c_test, E_test, m2_test = [], [], []
        grid_test = 'grid_beta_%s_test.txt' % beta
        with open(grid_test,'r') as g:
            for line in g:
                line = line.strip().split()
                c_test.append(line[0])
                E_test.append(int(line[1]))
                m2_test.append(float(line[2]))
        PlotEHist(E_list, E_bins, float(beta), test=E_test)
        PlotM2Hist(m2_list, m2_bins, float(beta), test=m2_test)
