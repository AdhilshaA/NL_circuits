import numpy as np
from math import exp
import matplotlib.pyplot as plt
from scipy.integrate import LSODA
figpath = r"H:\PERSONAL\PROJECTS\Programming\ACADEMICS\Int_lab2\NL_circuits\figures"



# Jerk Attractor



# Circuit parameters (resistances in Ohms, capacitances in Farads, voltages in Volts)
# FROM Msc thesis

R1 = 1e3
R2 = 1e3
R3 = 1e3
R4 = 1e3
R5 = 1e3
R6 = 1.167e3 #variable (potentiometer)
alpha = 0.026

C1 = 1e-6
C2 = 1e-6
C3 = 1e-6

# Initial conditions
state_0 = [0.0, 0.3, 0.0]
t_0 = 0.0
t_f = 0.5
max_step = 0.005

def equations(t, state):
    X, Y, Z = state
    
    # Non-linear circuit equations
    dxdt = Y / (C1 * R1)
    dydt = Z / (C2 * R2)
    dzdt = (- Z - ((R5/R6) * X) - (R5 * 1e-9 * (exp(Y / alpha) - 1))) / (C3 * R5)
    
    return [dxdt, dydt, dzdt]

# Solve the ODE system using LSODA method
sol = LSODA(equations,  t_0, state_0, t_f, vectorized=True)

# initialize the solution array
X = [ sol.y[0] ]
Y = [ sol.y[1] ]
Z = [ sol.y[2] ]
t = [ sol.t ]

print(sol.status)

while sol.status != 'finished':
    sol.step()
    t.append(sol.t)
    X.append(sol.y[0])
    Y.append(sol.y[1])
    Z.append(sol.y[2])
    

print(sol.status)

print(t[0], X[0], Y[0], Z[0])
print(t[-1], X[-1], Y[-1], Z[-1])
print(f"{len(t)} steps completed")

# # for creating a responsive plot
# %matplotlib widget
 
# creating figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(X, Y, Z, lw=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Jerk Attractor')
plt.show() 

# Plot XY, YZ, XZ projections of the attractor as 3 different plots and save them 
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(X, Y, lw=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('XY projection of Jerk Attractor')
plt.savefig(figpath+'\Jerk_sim_XY.png')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Y, Z, lw=0.5)
ax.set_xlabel('Z')
ax.set_ylabel('Y')
ax.set_title('YZ projection of Jerk Attractor')
plt.savefig(figpath+'\Jerk_sim_YZ.png')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(X, Z, lw=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_title('XZ projection of Jerk Attractor')
plt.savefig(figpath+'\Jerk_sim_XZ.png')

