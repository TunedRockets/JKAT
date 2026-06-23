import matplotlib.pyplot as plt
import numpy as np


ax = plt.figure().add_subplot(111, projection='3d', computed_zorder=True)

# Make data
r = 1
u = np.linspace(0, 2 * np.pi, 20)
v = np.linspace(0, np.pi, 20)
x = r * np.outer(np.cos(u), np.sin(v))
y = r * np.outer(np.sin(u), np.sin(v))
z = r * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='r', alpha=1)

r = np.linspace(0,2*np.pi,200)
x = np.cos(r)*3
y = np.sin(r)*3
z = r*0

ax.scatter(x,y,z)

plt.show()
