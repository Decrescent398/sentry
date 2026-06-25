import numpy as np

def angle_diff(a, b):
    return (a - b + np.pi) % (2*np.pi) - np.pi