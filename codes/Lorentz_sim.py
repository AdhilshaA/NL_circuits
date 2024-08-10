import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import LSODA

figpath = r"H:\PERSONAL\PROJECTS\Programming\ACADEMICS\Int_lab2\NL_circuits\figures"

# Lorenz Attractor


# Circuit parameters (resistances in Ohms, capacitances in Farads, voltages in Volts)
# FROM Msc thesis

R1 = 10e3
R2 = 100e3
R3 = 10e3
R4 = 100e3
R5 = 0.9765e3 #variable (potentiometer)
R6 = 5.6e3
R7 = 3.3e3
R8 = 3.6e3
R9 = 3.19e3
R10 = 100e3
R11 = 1e3
R12 = 3.3e3
R13 = 37.5e3
R14 = 3.3e3
R15 = 3.74e3
R16 = 100e3
R17 = 1e3
R18 = 1e3
R19 = 9e3

C1 = 200e-9
C2 = 200e-9
C3 = 200e-9

# Initial conditions
state_0 = [0.1, 0.1, 0.9]
t_0 = 0.0
t_f = 0.04
max_step = 0.005

def equations(t, state):
    X, Y, Z = state
    
    # Non-linear circuit equations
    dxdt = (- X - ((R4/R1) * X) + ((R4/R2) * X) + (R4/R3) * Y) / (C1 * R5)
    dydt = (- Y - ((R10/R7) *X * Z) + ((R10/R8) * X)) / (C2 * R11)
    dzdt = (- Z - ((R16/R13) * Z) + ((R16/R14) * X * Y)) / (C3 * R17)
    
    return [dxdt, dydt, dzdt]


# Solve the ODE system using LSODA method
sol = LSODA(equations,  t_0, state_0, t_f, vectorized=True)

# initialize the solution array
X = [sol.y[0]]
Y = [sol.y[1]]
Z = [sol.y[2]]
t = [sol.t]

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
ax.set_title('Lorenz Attractor')
plt.show()

# Plot XY, YZ, XZ projections of the attractor as 3 different plots and save them 
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(X, Y, lw=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('XY projection of Lorenz Attractor')
plt.savefig(figpath+'\Lor_sim_XY.png')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Z, Y, lw=0.5)
ax.set_xlabel('Z')
ax.set_ylabel('Y')
ax.set_title('YZ projection of Lorenz Attractor')
plt.savefig(figpath+'\Lor_sim_YZ.png')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(X, Z, lw=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_title('XZ projection of Lorenz Attractor')
plt.savefig(figpath+'\Lor_sim_XZ.png')

