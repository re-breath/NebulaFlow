import matplotlib.pyplot as plt
import numpy as np

# x = np.loadtxt('data.txt')[:,0]
# y = np.loadtxt('data.txt')[:,1]

x = [0, 1, 2, 3, 4, 5]
y = [0, 1, 4, 9, 16, 25]

#plt.plot(x, y)
plt.scatter(x, y, marker = 'h',color = 'r')

plt.show()

