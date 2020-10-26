import numpy as np
import pandas as pd
np.set_printoptions(precision=4)

def func(x, a):
    return a/x

pt = np.array([0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

names = [
    'dphiL123', 'dphiL124', 'dphiL134', 'dphiL234', 'dphiL12D1', 'dphiL13D1',
    'dphiL23D1', 'dphiL12D2', 'dphiL1D12', 'dphiL2D12', 'dphiL1D13', 'dphiL1D23',
    'dphiD123', 'dphiD124', 'dphiD134', 'dphiD234', 'dphiD235', 'dphiD245', 'dphiD345'
]

params = np.array([0.050, 0.073, 0.073, 0.050, 0.051, 0.070, 0.047, 0.059, 0.036, 0.036, 0.041, 0.041, 0.021, 0.031, 0.030, 0.020, 0.030, 0.030, 0.020])

for inum, iname in enumerate(names):
    print(iname)
    y = func(pt, params[inum])
    print(y, '\n')

