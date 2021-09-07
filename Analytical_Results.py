import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import json

mu = 1;
Y = []
Z = []
P = []

T = np.linspace(0, 1, 81);
X = np.linspace(0, 1, 33)

print(T)
print(X)

for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		y = (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5)
		Y.append(y)
	
with open('Y.txt', 'w') as f:
    for item in Y:
        f.write("%s\n" % item)
        
for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		z = (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5) * (1 + mu * math.pi**4) - 2 * mu * np.around(math.sin(x * math.pi), decimals=5)
		Z.append(z)

with open('Z.txt', 'w') as f:
    for item in Z:
        f.write("%s\n" % item)

for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		p = -1 * mu * ( (-1.0) * 2 * (1-t) * np.around(math.sin(x * math.pi), decimals=5) + math.pi**2 * (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5))
		P.append(p)

with open('P.txt', 'w') as f:
    for item in P:
        f.write("%s\n" % item)
