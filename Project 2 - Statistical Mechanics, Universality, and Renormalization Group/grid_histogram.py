#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:43:19 2019

@author: Minh Nguyen
"""

from matplotlib import pyplot as plt
from matplotlib import mlab
import numpy as np
import stats
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from math import log
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
    plt.bar(mid, n, width=width, align='center', label=r'data')
    #print(mid)
    E_mean, E_var, E_err, E_tau = stats.Stats(np.array(E))
    #print(E_ave, E_var, E_err, E_tau)
    plt.axvline(E_mean, lw=0.25, color='y', label=r'$\langle E \rangle$')
    
    #binWidth = (E_bin[-1] - E_bin[0])//len(E_bin)
    #x_arr = np.linspace(E_bin[0], E_bin[-1], 300, endpoint=True)
    #y_arr = Gaussian(x_arr, E_mean, E_var) * len(E)**(0.5)
    #x_arr += binWidth//2
    #plt.plot(x_arr, y_arr, color='r', label='curve')
    
    if len(test) > 0:
        PlotETest(test, E_bin, b)
    
    plt.errorbar(E_mean, max(n)/2, color='r', capsize=5, xerr=E_err, fmt='none', label=r'error')
    
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
    n,bins = np.histogram(M2, bins=100, density=True)
    #err = np.sqrt(n)#np.array([ np.sqrt(n[i]/E.count(bins[i])) for i in range(len(n)) ])
    width = (max(bins)-min(bins))/(len(bins)-1)
    mid = 0.5*(bins[1:] + bins[:-1])
    n *= width
    plt.bar(mid, n, width=width, align='center', label=r'data')
    M2_mean, M2_var, M2_err, M2_tau = stats.Stats(np.array(M2))
    #print(M2_mean, M2_var, M2_err, M2_tau)
    plt.axvline(M2_mean, lw=0.25, color='y', label=r'$\langle M^2\rangle$')
    
    #binWidth = (M2_bin[-1] - M2_bin[0])//len(M2_bin)
    #x_arr = np.array(sorted(list(set(abs(np.linspace(M2_bin[0], M2_bin[-1], 300, endpoint=True))))))
    #y_arr = Gaussian(x_arr, 0, M2_var) * len(M2)**(0.5)
    #x_arr += binWidth//2
    #plt.plot(x_arr, y_arr, color='r', label='curve')
    
    if len(test) > 0:
        PlotM2Test(test, bins, b)
        
    plt.errorbar(M2_mean, max(n)/2, color='r', capsize=5, xerr=M2_err, fmt='none', label=r'error')
    
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
    plt.bar(mid, n, width=width, align='center', label=r'test')
    E_mean, E_var, E_err, E_tau = stats.Stats(np.array(E))
    plt.axvline(E_mean, lw=0.25, label=r'test $\langle E \rangle$')

def PlotM2Test(M2, M2_bin, b):
    n,bins = np.histogram(M2, bins=M2_bin, density=True)
    width = 0.5*(max(bins)-min(bins))/(len(bins)-1)
    mid = 0.5*(bins[1:] + bins[:-1])
    n *= 2*width
    plt.bar(mid, n, width=width, align='center', label=r'test')
    M2_mean, M2_var, M2_err, M2_tau = stats.Stats(np.array(M2))
    plt.axvline(M2_mean, lw=0.25, label=r'test $\langle M^2\rangle$')

def PlotMatBeta(string, size, beta):
    if len(string) == size**2:
        str_split = [ s for s in string ]
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
    
def PlotMatCG(string, size, beta):
    if len(string) == size**2:
        str_split = [ s for s in string ]
        str_list = [ str_split[i:i+size] for i in range(0, len(str_split),size) ]
        str_list = np.array(str_list, int)
        fig = plt.figure()
        plt.title(r'Snapshot CG $%d\times%d$ at $\beta = %s$' % (size,size,beta))
        save_name = 'img/grid_snapshot_CG_%s_%d.png' % (beta, size)
        plt.imshow(str_list)
        plt.savefig(save_name, dpi=300)
    else:
        print('Cannot split string into a %dx%d grid for beta=%s' % (size, size, beta))

def PlotArrows(x0, y0, f, count=1):
    try:
        if count % 2 == 1:
            plt.arrow(x0, y0, f(x0)-x0, 0)
            PlotArrows(f(x0), y0, f, count=count+1)
        else:
            plt.arrow(x0, y0, 0, f(y0)-y0)
            PlotArrows(x0, f(y0), f, count=count+1)
    except ValueError:
        print('Arrows done')
    

'''    GRAPH E & M^2 BETA = 0.0, 0.1,..., 1.0, 1000000.0

beta_list = np.linspace(0,1,11, endpoint=True)
beta_list = [ '%.1f' % b for b in beta_list ]
beta_list.append('1000000.0')
#beta_list = [ str(b) for b in beta_list ]
#print(beta_list)
beta_test = ['0.0', '1000000.0']

m2_mean_list, m2_err_list = [], []

for beta in beta_list:
    c_list, E_list, m2_list = [], [], []
    #beta = beta_list[0]
    #print(beta)
    grid_file = 'grid_beta/grid_beta_%s.txt' % beta
    with open(grid_file,'r') as f:
        for line in f:
            line = line.strip().split()
            c_list.append(line[0])
            E_list.append(int(line[1]))
            m2_list.append(float(line[2]))
    
    
    #E_bins = np.arange(-2*len(c_list[0]), 2*len(c_list[0])+1, 4)
    #m2_bins = np.arange(0, max(m2_list)+0.01, 0.01)
    
    #PlotMat(c_list[-1], 27, float(beta))
    
    m2_mean, m2_var, m2_err, m2_tau = stats.Stats(np.array(m2_list))
    
    m2_mean_list.append(m2_mean)
    m2_err_list.append(m2_err)
    
    
    if beta in beta_test:
        c_test, E_test, m2_test = [], [], []
        grid_test = 'grid_beta/grid_beta_%s_test.txt' % beta
        with open(grid_test,'r') as g:
            for line in g:
                line = line.strip().split()
                c_test.append(line[0])
                E_test.append(int(line[1]))
                m2_test.append(float(line[2]))
        #PlotEHist(E_list, E_bins, float(beta), test=E_test, save=True)
        #PlotM2Hist(m2_list, m2_bins, float(beta), test=m2_test, save=True)
        E_mean,_,_,_ = stats.Stats(np.array(E_list))
        E_mean_test,_,_,_ = stats.Stats(np.array(E_test))
        print('beta =', beta,' dataE =', E_mean, ' testE =', E_mean_test)
        M2_mean,_,_,_ = stats.Stats(np.array(m2_list))
        M2_mean_test,_,_,_ = stats.Stats(np.array(m2_test))
        print('beta =', beta,' dataM2 =', M2_mean, ' testM2 =', M2_mean_test)
        print()
    else:
        #PlotEHist(E_list, E_bins, float(beta), save=True)
        #PlotM2Hist(m2_list, m2_bins, float(beta), save=True)
        E_mean,_,_,_ = stats.Stats(np.array(E_list))
        print('beta =', beta,' dataE =', E_mean)
        M2_mean,_,_,_ = stats.Stats(np.array(m2_list))
        print('beta =', beta,' dataM2 =', M2_mean)
        print()
'''



'''    PLOT M^2 vs beta

fig = plt.figure()
plt.title(r'$\langle M^2 \rangle$ vs $\beta$')
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\langle M^2 \rangle$')
plt.plot(beta_list[:-1], m2_mean_list[:-1], label='curve')
plt.errorbar(beta_list[:-1], m2_mean_list[:-1], color='r', capsize=5, 
             yerr=m2_err_list[:-1], fmt='none', label=r'error')
plt.legend()
plt.savefig('img/M2_vs_beta.png', dpi=300)
'''


'''    PLOT C_v vs T and E vs T

E_mean_list, T_list = [], []

Cv_file = 'grid_Cv.txt'
with open(Cv_file, 'r') as f:
    for line in f:
        line = np.array(line.strip().split(), dtype=float)
        E_mean, E_var, E_err, E_tau = stats.Stats(line[:-1])
        
        E_mean_list.append(E_mean)
        T_list.append(line[-1])

#b_list = np.array(T_list)**(-1)

fig = plt.figure()
plt.title(r'$\langle E \rangle$ vs $T$')
plt.xlabel(r'$T$')
plt.ylabel(r'$\langle E \rangle$ and $C_v$')
plt.plot(T_list, E_mean_list, label=r'$\langle E \rangle$ vs $T$')
#plt.plot(b_list, E_mean_list, label=r'$\langle E \rangle$ vs $\beta$')

T_list.sort()
#b_list.sort()
E_mean_list.reverse()

spl_T = UnivariateSpline(T_list, E_mean_list).derivative()

Cv = list(spl_T(T_list))
plt.plot(T_list, Cv, label=r'$C_v = \frac{\partial E}{\partial T}$ vs $T$')

ind_max_Cv = Cv.index(max(Cv))
plt.axvline(T_list[ind_max_Cv], color='y', label=r'Max $C_v$ at $T = %.6f$' % T_list[ind_max_Cv])

plt.legend()
plt.savefig('img/grid_Cv_vs_T.png', dpi=300)
'''

'''    PLOT TEST CG SNAPSHOTS
file_name = 'grid_rand_CG_compare.txt'
with open(file_name, 'r') as f:
    for line in f:
        line = line.strip().split()
        state = line[0]
        size = int(line[1])
        beta = '0.0'
        
        str_split = [ s for s in state ]
        str_list = [ str_split[i:i+size] for i in range(0, len(str_split),size) ]
        str_list = np.array(str_list, int)
        fig = plt.figure()
        plt.title(r'Snapshot test $%d\times%d$ at $\beta = %s$' % (size,size,beta))
        save_name = 'grid_snapshot_CG_%s_%d_test.png' % (beta, size)
        plt.imshow(str_list)
        plt.savefig(save_name, dpi=300)
'''

'''    Plot coarse grained snapshots comparison 

beta_list, state_list = [], []
file_name = 'grid_RG/grid_RG_compare_cg.txt'
with open(file_name, 'r') as f:
    for line in f:
        line = line.strip().split()
        beta = '%.1f' % float(line[0])
        state_list.append([beta, line[1:]])

for state in state_list:
    beta = state[0]
    for s in state[1]:
        size = int(round(len(s)**0.5))
        PlotMatCG(s, size, beta)
'''


'''    PLOT CG VS NATIVE CURVES FOR M^2 VS beta    

nat_m2_mean_list, cg_m2_mean_list = [], []
beta_list = []
test_m2_nat, test_m2_cg = [], []

native_file = 'grid_RG/grid_RG.txt'
cg_file     = 'grid_RG/grid_RG_cg.txt'
test_file   = 'grid_RG/grid_RG_0.0_compare_cg_test.txt'

with open(native_file,'r') as nat, open(cg_file,'r') as cg:
    for natline in nat:
        natline = natline.strip().split()
        beta_list.append(float(natline[0]))
        
        nat_m2 = natline[1:]
        nat_m2 = np.array(nat_m2, float)
        nat_m2_mean,_,_,_ = stats.Stats(nat_m2)
        nat_m2_mean_list.append(nat_m2_mean)
    
    for cgline in cg:
        cgline = cgline.strip().split()
        cg_m2 = cgline[1:]
        cg_m2 = np.array(cg_m2, float)
        cg_m2_mean,_,_,_ = stats.Stats(cg_m2)
        cg_m2_mean_list.append(cg_m2_mean)
        
with open(test_file,'r') as test:
    for line in test:
        line = line.strip().split()
        test_m2_nat.append(float(line[0]))
        test_m2_cg.append(float(line[1]))

test_m2_nat = np.array(test_m2_nat)
test_m2_cg  = np.array(test_m2_cg)

test_mean_m2_nat,_,_,_ = stats.Stats(test_m2_nat)
test_mean_m2_cg,_,_,_ = stats.Stats(test_m2_cg)

fig = plt.figure()
plt.title(r'$\langle M^2 \rangle$ vs $\beta$')
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\langle M^2 \rangle$')
plt.plot(beta_list, nat_m2_mean_list, label=r'Native (nat)')
plt.plot(beta_list, cg_m2_mean_list, label=r'Coarse grained (CG)')
plt.plot(0.0, test_mean_m2_nat,'o', 
         label=r'Test (nat) $\langle M^2 \rangle= %f$' % test_mean_m2_nat)
plt.plot(0.0, test_mean_m2_cg,'o', 
         label=r'Test (CG) $\langle M^2 \rangle= %f$' % test_mean_m2_cg)
plt.legend()
plt.savefig('M2_vs_beta_native_vs_CG.png', dpi=300)
#plt.show()
'''
    

'''    PLOTTING R(J) vs J    '''
nat_m2_mean_list, cg_m2_mean_list = [], []
beta_list = []

native_file = 'grid_RG/grid_RG.txt'
cg_file     = 'grid_RG/grid_RG_cg.txt'

with open(native_file,'r') as nat, open(cg_file,'r') as cg:
    for natline in nat:
        natline = natline.strip().split()
        beta_list.append(float(natline[0]))
        
        nat_m2 = natline[1:]
        nat_m2 = np.array(nat_m2, float)
        nat_m2_mean,_,_,_ = stats.Stats(nat_m2)
        nat_m2_mean_list.append(nat_m2_mean)
    
    for cgline in cg:
        cgline = cgline.strip().split()
        cg_m2 = cgline[1:]
        cg_m2 = np.array(cg_m2, float)
        cg_m2_mean,_,_,_ = stats.Stats(cg_m2)
        cg_m2_mean_list.append(cg_m2_mean)

arr = np.linspace(0.01,0.998,1000)

nat_R_func   = interp1d(nat_m2_mean_list, beta_list, fill_value='extrapolate')

nat_R_J     = list(nat_R_func(arr))
R_J         = nat_R_func(cg_m2_mean_list)

fig = plt.figure()
plt.title(r'$R(J)$ vs $J$')
plt.xlabel(r'$J$')
plt.ylabel(r'$R(J)$')
plt.plot(beta_list, R_J, label=r'$R(J)$ vs $J$')
x = np.linspace(0, 1.4, 1000)
y = x
plt.plot(x,y, '--', label=r'y=x')


# Calculate critical point
cg_R_func = interp1d(cg_m2_mean_list, beta_list)
cg_R_J = list(cg_R_func(arr))
crit_trans_T = 0
temp = 1
for nat, cg in zip(nat_R_J, cg_R_J):
    if abs(nat-cg) < temp:
        crit_trans_T = cg
        temp = abs(nat-cg)

plt.plot(crit_trans_T, crit_trans_T, 'o', 
         label=r'Critical beta $\approx$ %.3f' % crit_trans_T)


# Graphing arrows and calculate critical exponent
R = interp1d(beta_list, R_J)#, fill_value='extrapolate')
dR_dJ = UnivariateSpline(beta_list, R_J).derivative()

#PlotArrows(0.35, R(0.35), R)
#PlotArrows(0.5, R(0.5), R)

crit_trans_T = beta_list[list(dR_dJ(beta_list)).index(max(dR_dJ(beta_list)))]
print(crit_trans_T)
crit_exp = log(dR_dJ(crit_trans_T),3.0)
print(crit_exp)

plt.legend()
#plt.show()
plt.savefig('img/RJ_vs_J_no_arrows.png', dpi=300)

















