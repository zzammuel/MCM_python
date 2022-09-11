import numpy as np
import pandas as pd

from math import lgamma

####################
### Log-evidence ###
####################

def complexity(m):
    power = 1 << (m-1)
    return np.log(np.pi) * power - np.math.lgamma(power)

def project_icc(Kset, ICC):
    new_Kset = {}
    for s,ks in Kset.items():
        s_new = s & ICC
        if s_new in new_Kset:
            new_Kset[s_new] += ks
        else:
            new_Kset[s_new] = ks
    return new_Kset, bin(ICC).count("1")

def log_evidence_ICC(Kset, ICC, N):
    LogE = 0
    new_Kset, nNode = project_icc(Kset, ICC)
    for s,ks in new_Kset.items():
        LogE += np.math.lgamma(ks + 0.5)
    LogE -= (complexity(nNode) + np.math.lgamma(N + (1 << (nNode - 1))))
    return LogE

def log_evidence_MCM(Kset, MCM, N):
    LogE = 0
    for _,ICC in MCM.items():
        LogE += log_evidence_ICC(Kset, ICC, N)
    return LogE

######################
### Log-likelihood ###
######################

def log_likelihood_ICC(Kset, ICC, N):
    LogL = 0
    new_Kset, nNode = project_icc(Kset, ICC)
    for s,ks in new_Kset.items():
        LogL += ks * np.math.log(ks/N)
    return LogL

def log_likelihood_MCM(Kset, MCM, N):
    LogL = 0
    for _,ICC in MCM.items():
        LogL += log_likelihood_ICC(Kset, ICC, N)
    return LogL

#########################
### Merging procedure ###
#########################

def merging(Kset, N, n):
    MCM = {}
    logE_MCM = {}
    
    for i in range(n):
        MCM[i] = (1 << i)
        logE_MCM[i] = log_evidence_ICC(Kset, MCM[i], N)
    
    logE_mat = np.zeros((n,n))
    
    stop = False
    
    while not stop and len(MCM)-1:
        best_diff = 0
        for i,node1 in MCM.items():
            for j,node2 in MCM.items():
                if i>=j:
                    continue

                if not logE_mat[i,j]:
                    logE_mat[i,j] = log_evidence_ICC(Kset, node1 + node2, N)
                    logE_mat[j,i] = logE_mat[i,j]
                
                if logE_mat[i,j] - logE_MCM[i] - logE_MCM[j] > best_diff:
                    best_diff = logE_mat[i,j] - logE_MCM[i] - logE_MCM[j]
                    best_merger = [[i,node1], [j,node2]]
                    
        if best_diff == 0:
            stop = True
            continue
        
        i_keep = best_merger[0][0]; node_keep = best_merger[0][1]
        i_erase = best_merger[1][0]; node_erase = best_merger[1][1]
                
        logE_MCM[i_keep] = logE_mat[i_keep, i_erase] 
        MCM[i_keep] = node_keep + node_erase
        
        logE_MCM.pop(i_erase)
        MCM.pop(i_erase)
        
        logE_mat[i_keep,:] = 0
        logE_mat[:,i_keep] = 0
        
    return MCM, sum([items for _,items in logE_MCM.items()])