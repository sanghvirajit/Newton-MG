import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

y = pd.read_csv("Residual.txt")
x = pd.read_csv("Iterations.txt")

fig = plt.gcf()

plt.semilogy(x, y, label="Residual")
plt.grid()
plt.xlabel("SOR Iterations")
plt.ylabel("Residual")
fig.savefig("Residual.png")
plt.show()


