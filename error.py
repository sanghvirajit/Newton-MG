import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from numpy import linalg as LA

x = 0.9375

y_df = pd.read_csv("Y.txt", header=None)

yh = pd.read_csv("y_experimental_{}.txt".format(x), header=None)
y = pd.DataFrame(index=np.arange((yh.size)), columns=np.arange(1))

n = 14

for i in range(0, 257):
	y.iloc[i, 0] = y_df.iloc[i*15 + n, 0]

y_diff = y - yh

y_array = y_diff.to_numpy()

print("y error: {}".format(LA.norm(y_array)))

p_df = pd.read_csv("P.txt", header=None)

ph = pd.read_csv("p_experimental_{}.txt".format(x), header=None)
p = pd.DataFrame(index=np.arange((ph.size)), columns=np.arange(1))

for i in range(0, 257):
	p.iloc[i, 0] = p_df.iloc[i*15 + n, 0]

p_diff = p - ph

p_array = p_diff.to_numpy()

print("p error: {}".format(LA.norm(p_array)))
