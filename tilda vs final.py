import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

mu = 1;
x = 0.03125;

t = np.linspace(0, 1, 81)

#PLOT Y
y1 = (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5)

y2 = pd.read_csv("y_tilda_experimental_{}.txt".format(x), header=None)
y3 = pd.read_csv("y_tilda_vector_experimental_{}.txt".format(x), header=None)
y4 = pd.read_csv("y_final_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()
plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

#plt.plot(t, y3, "k+")
plt.plot(t, y2, label="y tilda")
plt.plot(t, y3, label="y tilda vector")
plt.plot(t, y4, label="y final")

#plt.plot(t, y2, "kx")
plt.plot(t, y1, label="y analytical")

plt.xlabel("Time step")
plt.ylabel("y")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical Y_{}.png".format(x))
plt.show()

#PLOT P
p1 = ( 2 * mu * (1-t) * np.around(math.sin(x * math.pi), decimals=5) - mu * math.pi**2 * (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5))
p2 = pd.read_csv("p_tilda_experimental_{}.txt".format(x), header=None)
p3 = pd.read_csv("p_tilda_vector_experimental_{}.txt".format(x), header=None)
p4 = pd.read_csv("p_final_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

#plt.plot(t, p1, "k+")
plt.plot(t, p2, label="p tilda")
plt.plot(t, p3, label="p tilda vector")
plt.plot(t, p4, label="p final")

#plt.plot(t, p2, "kx")
plt.plot(t, p1, label="p analytical")

plt.xlabel("Time step")
plt.ylabel("p")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical P_{}.png".format(x))
plt.show()


