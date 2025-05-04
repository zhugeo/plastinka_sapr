import csv
from matplotlib import pyplot as plt

with open("build/innerNodes.csv") as csvfile:
    raw = csv.reader(csvfile, delimiter=" ", quotechar="|")
    innerNodes = list(raw)

with open("build/outerNodes.csv") as csvfile:
    raw = csv.reader(csvfile, delimiter=" ", quotechar="|")
    outerNodes = list(raw)

X_inner = [float(i[0]) for i in innerNodes[1:]]
Y_inner = [float(i[1]) for i in innerNodes[1:]]

X_outer = [float(i[0]) for i in outerNodes[1:]]
Y_outer = [float(i[1]) for i in outerNodes[1:]]

plt.plot(X_inner, Y_inner, "bo")
plt.plot(X_outer, Y_outer, "r+")

plt.show()
