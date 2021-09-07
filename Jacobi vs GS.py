import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

y1 = pd.read_csv("Residual_Jacobi.txt")
x1 = pd.read_csv("Iterations_Jacobi.txt")

y2 = pd.read_csv("Residual_GaussSeidel.txt")
x2 = pd.read_csv("Iterations_GaussSeidel.txt")

y3 = pd.read_csv("Residual_SOR.txt")
x3 = pd.read_csv("Iterations_SOR.txt")

fig = plt.gcf()

plt.semilogy(x1, y1, label="Jacobi")
plt.semilogy(x2, y2, label="Gauss Seidel")
plt.semilogy(x3, y3, label="SOR")
plt.grid()

plt.xlabel("Iterations")
plt.ylabel("Residual")
plt.legend()

fig.savefig("Residual Jacobi Vs GS Vs SOR.png")
plt.show()
