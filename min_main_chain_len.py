from math import log,ceil

def min_main_chain_len(NUM_ATOMS):
    return min(2*ceil(log((NUM_ATOMS-1)/2+1,3))+1, 2*ceil(log(NUM_ATOMS+1,3)))

from sys import argv

print(min_main_chain_len(int(argv[1])))