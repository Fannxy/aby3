import numpy as np
import math

def get_k(V, E, B=1024):
    B = max(B, V)
    k = int(np.ceil(B*V/E))
    k = 2**math.floor(math.log2(k))
    
    return k
