import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

T = np.linspace(0, 1, 81)

u = pd.read_csv("Z.txt", header=None)
df = pd.DataFrame(index=np.arange((u.size)/31), columns=np.arange(1))

n = 0

for i in range(0, 81):
	df.iloc[i, 0] = u.iloc[i*31 + n, 0]
	
df.to_csv("Z_final.txt", header=None, index=False)

fig = plt.gcf()

plt.plot(T, df, label="z")

plt.show()
