## Read a binary datafile
## Data file should be strings of 0s and 1s without spacing
def read_data(filename, n, flip = False):
    Kset = {}
    N = 0
    with open(filename, 'r') as file:
        for line in file.readlines():
            s = int(line[:n],2)
            if flip:
                s ^= 2**n - 1
            N += 1
            if s in Kset:
                Kset[s] += 1
            else:
                Kset[s] = 1
    return Kset, N

## Basis transformation
def transform(Nset, basis):
    Kset = {}
    
    for s, ks in Nset.items():
        Op = 1
        sn = 0
        for nr, op in basis.items():
            if not bin(s & op).count("1")%2:
                sn += Op
            Op <<= 1
        if sn in Kset:
            Kset[sn] += ks
        else:
            Kset[sn] = ks

    return Kset

## Transform a partition to the c1_c_2_..._c_n format
def rewrite_partition(MCM, n):
    string = ""
    for i in range(n):
        for key, val in MCM.items():
            binary = bin(val).replace("0b", "")
            binary = '0' * (n - len(binary)) + binary
            if binary[i] == '1':
                string += str(key) + "_"
                break
    return string[:-1]

## Transform a partition to a dict of binary representations of the communities
def restore_partition(MCM, n):
    MCM_rewritten = {}
    Op = 1 << n - 1

    for label in MCM.split('_'):
        if label in MCM_rewritten:
            MCM_rewritten[label] += Op
        else:
            MCM_rewritten[label] = Op
        Op >>= 1
    return MCM_rewritten

## Print a partition in terms of its binary representation
def print_partition(MCM, n):
    for com,ICC in MCM.items():
        binary = bin(ICC).replace("0b", "")
        missing = n - len(binary)
        print(f'com {com}\t' + '0'*missing + binary)