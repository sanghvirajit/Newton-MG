import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

mu = 0.1;
x = 0.0625;

t = np.linspace(0, 1, 257)

#PLOT Y
y1 = pd.read_csv("y_final_experimental_{}.txt".format(x), header=None)
y2 = (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5)

fig = plt.gcf()
plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

#plt.plot(t, y1, "k+")
plt.plot(t, y1, label="y exp")

#plt.plot(t, y2, "kx")
plt.plot(t, y2, label="y analytical")

plt.xlabel("Time step")
plt.ylabel("y")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical Y_{}.png".format(x))
plt.show()

#PLOT P
p1 = pd.read_csv("p_final_experimental_{}.txt".format(x), header=None)
p2 = ( 2 * mu * (1-t) * np.around(math.sin(x * math.pi), decimals=5) - mu * math.pi**2 * (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5))

fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

#plt.plot(t, p1, "k+")
plt.plot(t, p1, label="p exp")

#plt.plot(t, p2, "kx")
plt.plot(t, p2, label="p analytical")

plt.xlabel("Time step")
plt.ylabel("p")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical P_{}.png".format(x))
plt.show()


