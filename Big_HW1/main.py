import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Function to find y(x)
def findy(x):
    return np.sin(3 * x + 0.2)

# Function for finding length of the plot
def integrate(x,a,b):
    return np.sqrt(1+(a**2)*(np.cos(a*x+b))**2)  # length of sin(3x+0.2) is integral of
                                                 # sqrt(1+((sin(3x+0.2))')^2) = sqrt(1+9(cos(3x+0.2))^2)

# Find radius of curvacy at given point
def findR(x):
    return np.abs(np.sqrt((1 + (3 * np.cos(3 * x + 0.2)) ** 2) ** 3) / (-9 * np.sin(3 * x + 0.2)))

# Euclidean distance
def find_distance(x1, x2):
    y1 = findy(x1)
    y2 = findy(x2)
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)

# Function to find normal acceleration
def find_an(x, v):
    R = findR(x)
    return (v ** 2) / R

# Required constraints
dx = 0.0001
da = 0.1
x_data = np.arange(0, 4 + dx / 2, dx)
y_data = []
v_data = np.zeros(len(x_data) + 1, dtype=float)
v_max = 1.5
a_t_max = 10
a_t_real = np.full(len(x_data), a_t_max, dtype=float)
a_n_max = 6
a_n_real = np.full(len(x_data), 0, dtype=float)
d_data = np.zeros(len(x_data), dtype=float)

for i in range(len(x_data)):
    y_data.append(findy(x_data[i]))

i = 0
ind = np.where(x_data == 3.5)
index = ind[0][0]

# Find velocity for 0 <= x < 3.5
while(i < index):
    distance = find_distance(x_data[i], x_data[i + 1])
    d_data[i] = distance
    a_t_cur = a_t_real[i]
    v_cur = v_data[i]
    v_fut = np.sqrt(2 * distance * a_t_cur + v_cur ** 2)
    a_n_cur = find_an(x_data[i + 1], v_fut)
    if (v_fut > v_max or a_n_cur > a_n_max):
        if (a_t_real[i] - da < -a_t_max):
            a_t_real[i] = a_t_max
            a_t_real[i - 1] -= da
            i -= 1
        else:
            a_t_real[i] -= da
        continue
    else:
        a_t_real[i] = a_t_cur
        a_n_real[i] = a_n_cur
        v_data[i + 1] = v_fut
    i += 1

# Find velocity for 3.5 <= x < 4
j = len(x_data) - 1
while(j >= index):
    distance = find_distance(x_data[j], x_data[j - 1])
    d_data[j] = distance
    a_t_cur = a_t_real[j]
    v_cur = v_data[j]
    v_fut = np.sqrt(2 * distance * a_t_cur + v_cur ** 2)
    a_n_cur = find_an(x_data[j - 1], v_fut)
    if (v_fut > v_max or a_n_cur > a_n_max):
        if (a_t_real[j] - da < -a_t_max):
            a_t_real[j] = a_t_max
            a_t_real[j + 1] -= da
            j += 1
        else:
            a_t_real[j] -= da
        continue
    else:
        a_t_real[j] = a_t_cur
        a_n_real[j] = a_n_cur
        v_data[j - 1] = v_fut
    j -= 1

for i in range(index, len(a_t_real)):
    a_t_real[i] = -a_t_real[i]

# Calculate the total time
t = 0
v_data = np.delete(v_data, 0)
for i in range(len(v_data) - 2):
    t_cur = d_data[i] / v_data[i]
    t += t_cur

print("Total time is:", t, "\n")
print("The whole distance as sum of short sectors:", np.sum(d_data))
print("The whole distance using integral of function:", quad(integrate, 0, 4, args=(3, 0.2))[0])

# Plot y(x)
plt.plot(x_data, y_data)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# Plot y(t)
plt.plot(np.arange(0, t, t / len(y_data)), y_data)
plt.xlabel('Time')
plt.ylabel('y')
plt.show()

# Plot v(t)
plt.plot(np.arange(0, t, t / len(y_data)), v_data)
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.show()

# Plot tangential acceleration
plt.plot(np.arange(0, t, t / len(a_t_real)), a_t_real)
plt.xlabel('Time')
plt.ylabel('Tangential acceleration')
plt.show()

# Plot normal acceleration
plt.plot(np.arange(0, t, t / len(a_t_real)), a_n_real)
plt.xlabel('Time')
plt.ylabel('Normal acceleration')
plt.show()