import numpy as np
import matplotlib.pyplot as plt

# Load the data file
data = np.loadtxt('outputs/test.txt')

# Extract columns
x1 = data[:, 0]
y1 = data[:, 1]
x2 = data[:, 2]
y2 = data[:, 3]

# Create the figure and plot segments
plt.figure(figsize=(8, 8))

# First loop: Plot blue lines between (x1, y1) and (x2, y2)
for i in range(len(x1)):
    plt.plot([x1[i], x2[i]], [y1[i], y2[i]], color='blue' , linestyle='-', marker='o', 
             markersize=2, markerfacecolor='black', markeredgecolor='blue', linewidth=1.5)

# Second loop: Plot red normal vector arrows for a subset of the data
limit = int(len(x1) / 10)
for j in range(limit):
    i = 10 * j
    halfx = (x2[i] - x1[i]) / 2
    halfy = (y2[i] - y1[i]) / 2
    norm_normale = ((y2[i] - y1[i])**2 + (x1[i] - x2[i])**2)**(1/2)
    
    # Avoid division by zero if segment length is zero
    if norm_normale == 0:
        continue

    # Start point (midpoint of the segment)
    start_x = x1[i] + halfx
    start_y = y1[i] + halfy
    
    # Direction components (u, v) for the normal vector
    u = (y2[i] - y1[i]) / norm_normale
    v = (x1[i] - x2[i]) / norm_normale

    # Use quiver to draw the arrow
    plt.quiver(start_x, start_y, u, v, angles='xy', scale_units='xy', scale=1, color='red', width=0.005)

# Formatting
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Line segments with normal vector arrows")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')

# Force matplotlib to update data limits to include everything and add a 10% margin padding
plt.autoscale()
x_root, x_high = plt.xlim()
y_root, y_high = plt.ylim()
x_pad = (x_high - x_root) * 1
y_pad = (y_high - y_root) * 1
plt.xlim(x_root - x_pad, x_high + x_pad)
plt.ylim(y_root - y_pad, y_high + y_pad)

# Display the plot
plt.show()