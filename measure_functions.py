import numpy as np

def norm_mut_info(MCM1, MCM2, n):
    I = 0
    H = 0
    flag = 0
    for _,ICC1 in MCM1.items():
        p1 = bin(ICC1).count("1") / n
        for _,ICC2 in MCM2.items():
            p2 = bin(ICC2).count("1") / n
            p12 = bin(ICC1 & ICC2).count("1") / n
            if p12:
                I += p12 * np.log(p12/(p1*p2))
            if flag < len(MCM2):
                H += p2 * np.log(p2)
                flag += 1
        H += p1 * np.log(p1)
    
    if not H:
        return 1
    else:
        return -2*I/H

def var_of_info(MCM1, MCM2, n):
    I = 0
    H = 0
    for _,ICC1 in MCM1.items():
        p1 = bin(ICC1).count("1") / n
        for _,ICC2 in MCM2.items():
            p2 = bin(ICC2).count("1") / n
            p12 = bin(ICC1 & ICC2).count("1") / n
            if p12:
                I += p12 * np.log(p12/(p1*p2))
                H += p12 * np.log(p12)
    
    if not H:
        return 0
    else:
        return 1 + I/H