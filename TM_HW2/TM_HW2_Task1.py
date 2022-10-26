import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Create plot
fig, ax = plt.subplots()
ax.set_xlim(-100, 120)
ax.set_ylim(-100, 80)
ax.set_aspect("equal")
ax.grid()

# Create lines for the whole system and for each vector of velocity and acceleration
line, = ax.plot([], [], 'o-', lw=2)

line_vA, = ax.plot([], [], '<-', lw=1, color="blue", label="vA")
line_aA, = ax.plot([], [], '<-', lw=1, color="black", label="aA")

line_vB, = ax.plot([], [], '<-', lw=1, color="green", label="vB")
line_aB, = ax.plot([], [], '<-', lw=1, color="red", label="aB")

line_vC, = ax.plot([], [], '<-', lw=1, color="orange", label="vC")

line_vD, = ax.plot([], [], '<-', lw=1, color="pink", label="vD")

line_vE, = ax.plot([], [], '<-', lw=1, color="purple", label="vE")

line_vF, = ax.plot([], [], '<-', lw=1, color="yellow", label="vF")

# Text where angular velocity will be shown
text = ax.text(0.2, 0.25, 'matplotlib', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
leg = ax.legend(loc='upper left')

# Initial given values
w = 2
phi0 = np.pi/3
AO1 = 21
BO2 = 25
FO3 = 20
AB = 54
BC = 52
CD = 69
CE = 35
EF = 32
AC = CD/3

# Coordinates of O1, O2 and O3
xO1, yO1 = 0, 0
xO2, yO2 = 56, -26
xO3, yO3 = 66, 41

# Draw circles that are trajectory of points A, B and F
ax.add_artist(plt.Circle((xO1, yO1), AO1, linestyle='--', color="None", ec="orange"))
ax.add_artist(plt.Circle((xO2, yO2), BO2, linestyle='--', color="None", ec="blue"))
ax.add_artist(plt.Circle((xO3, yO3), FO3, linestyle='--', color="None", ec="red"))


# Arrays for coordinates of each point
Ax, Ay, vAx, vAy, aAx, aAy = [], [], [], [], [], []
Bx, By, vBx, vBy, aBx, aBy = [], [], [], [], [], []
Cx, Cy, vCx, vCy = [], [], [], []
Dx, Dy, vDx, vDy = [], [], [], []
Ex, Ey, vEx, vEy = [], [], [], []
Fx, Fy, vFx, vFy = [], [], [], []

# All required angles
phi, beta, alpha, sigma, dzeta, ksi, eta, gamma = [], [], [], [], [], [], [], []

# Time from 0 to 2pi with step 0.01
dt = 0.01
t = np.arange(0, np.pi, dt)

# Necessary constants
O1O2 = np.sqrt((xO1-xO2)**2+(yO1-yO2)**2)
theta = np.arccos(xO2/O1O2)
ro = np.arccos((AC**2 + AB**2 - BC**2)/(2*AC*AB))

# Functions to find coordinates of each point
# Using derived angles and coordinates of other points
# we can find coordinates of each point with respect to time
def findA():
    for i in range(len(t)):
        phi_cur = phi0 + w * t[i]
        phi.append(phi_cur)
        Ax.append(AO1 * np.cos(phi_cur))
        Ay.append(AO1 * np.sin(phi_cur))

def findB():
    for i in range(len(t)):
        AO2 = np.sqrt((xO2-Ax[i])**2+(yO2-Ay[i])**2)

        beta_cur = np.arcsin((AO1/AO2)*np.sin(phi[i]+theta))
        beta.append(beta_cur)

        alpha_cur = np.arccos((BO2**2 + AO2**2 - AB**2)/(2*BO2*AO2))
        alpha.append(alpha_cur)

        angle = np.pi - (alpha_cur + beta_cur + theta)

        Bx.append(xO2 + BO2 * np.cos(angle))
        By.append(yO2 + BO2 * np.sin(angle))

def findC():
    for i in range(len(t)):
        sigma1 = np.pi - (phi[i] + theta + beta[i])
        sigma2 = np.arcsin(BO2/AB * np.sin(alpha[i]))
        sigma_cur = sigma1 + sigma2
        sigma.append(sigma_cur)

        dz_cur = np.pi - (sigma_cur + phi[i])
        dzeta.append(dz_cur)

        ksi_cur = ro - dz_cur
        ksi.append(ksi_cur)

        Cx.append(Ax[i] + AC * np.cos(ksi_cur))
        Cy.append(Ay[i] + AC * np.sin(ksi_cur))

def findD():
    for i in range(len(t)):
        Dy.append(16)

        eta_cur = np.arcsin((Cy[i]-Dy[i])/CD)
        eta.append(eta_cur)
        Dx.append(Cx[i] + CD * np.cos(eta_cur))

def findE():
    for i in range(len(t)):
        Ex.append(Cx[i] + CE * np.cos(eta[i]))
        Ey.append(Cy[i] - CE * np.sin(eta[i]))

def findF():
    for i in range(len(t)):
        EO3 = np.sqrt((Ex[i] - xO3)**2 + (Ey[i] - yO3)**2)
        CO3 = np.sqrt((Cx[i] - xO3) ** 2 + (Cy[i] - yO3) ** 2)

        eta1 = np.arccos((CE**2 + EO3**2 - CO3**2)/(2*CE*EO3))
        eta2 = np.arccos((EF**2 + EO3**2 - FO3**2)/(2*EF*EO3))

        gamma_cur = np.pi - (eta1 + eta2 + eta[i])
        gamma.append(gamma_cur)

        Fx.append(Ex[i] + EF * np.cos(gamma_cur))
        Fy.append(Ey[i] + EF * np.sin(gamma_cur))

# Function to find angular velocity of link using fomula: |w_AB| = |v_AB| / AB
def findw(vx1, vy1, vx2, vy2, l):
    w = []
    for i in range(len(t)):
        w.append(np.sqrt((vx1[i]-vx2[i])**2 + (vy1[i]-vy2[i])**2)/l)
    return w

# Find all coordinates
findA()
findB()
findC()
findD()
findE()
findF()

# Functions to find velocities and accelerations of each point using derivatives
def find_vA_aA():
    global vAx, vAy, aAx, aAy
    vAx = np.diff(Ax) / dt
    vAy = np.diff(Ay) / dt
    vAx = np.insert(vAx, 0, 0)
    vAy = np.insert(vAy, 0, 0)
    aAx = np.diff(vAx) / dt
    aAy = np.diff(vAy) / dt
    aAx = np.append(aAx, 0)
    aAy = np.append(aAy, 0)

def find_vB_aB():
    global vBx, vBy, aBx, aBy
    vBx = np.diff(Bx) / dt
    vBy = np.diff(By) / dt
    vBx = np.insert(vBx, 0, 0)
    vBy = np.insert(vBy, 0, 0)
    aBx = np.diff(vBx) / dt
    aBy = np.diff(vBy) / dt
    aBx = np.append(aBx, 0)
    aBy = np.append(aBy, 0)

def find_vC():
    global vCx, vCy
    vCx = np.diff(Cx) / dt
    vCy = np.diff(Cy) / dt
    vCx = np.insert(vCx, 0, 0)
    vCy = np.insert(vCy, 0, 0)

def find_vD():
    global vDx, vDy
    vDx = np.diff(Dx) / dt
    vDy = np.diff(Dy) / dt
    vDx = np.insert(vDx, 0, 0)
    vDy = np.insert(vDy, 0, 0)

def find_vE():
    global vEx, vEy
    vEx = np.diff(Ex) / dt
    vEy = np.diff(Ey) / dt
    vEx = np.insert(vEx, 0, 0)
    vEy = np.insert(vEy, 0, 0)

def find_vF():
    global vFx, vFy
    vFx = np.diff(Fx) / dt
    vFy = np.diff(Fy) / dt
    vFx = np.insert(vFx, 0, 0)
    vFy = np.insert(vFy, 0, 0)

# Find velocities and accelerations of each point
find_vA_aA()
find_vB_aB()
find_vC()
find_vD()
find_vE()
find_vF()


# Find angular velocity of each link
wAO1 = findw(vAx, vAy, np.zeros(len(t)), np.zeros(len(t)), AO1)
wBO2 = findw(vBx, vBy, np.full(len(t), xO2), np.full(len(t), yO2), BO2)
wAB = findw(vAx, vAy, vBx, vBy, AB)
wAC = findw(vAx, vAy, vCx, vCy, AC)
wBC = findw(vBx, vBy, vCx, vCy, BC)
wCE = findw(vEx, vEy, vCx, vCy, CE)
wCD = findw(vDx, vDy, vCx, vCy, CD)
wEF = findw(vEx, vEy, vFx, vFy, EF)
wFO3 = findw(vFx, vFy, np.full(len(t), xO3), np.full(len(t), yO3), FO3)

# Animation function
def animation_frame(i):
    line.set_data([xO1, Ax[i], Bx[i], xO2, Bx[i], Cx[i], Ax[i], Cx[i], Dx[i], Ex[i], Fx[i], xO3],
                  [yO1, Ay[i], By[i], yO2, By[i], Cy[i], Ay[i], Cy[i], Dy[i], Ey[i], Fy[i], yO3])

    line_vA.set_data([Ax[i], Ax[i] + vAx[i]], [Ay[i], Ay[i] + vAy[i]])
    line_aA.set_data([Ax[i], Ax[i] + aAx[i]], [Ay[i], Ay[i] + aAy[i]])

    line_vB.set_data([Bx[i], Bx[i] + vBx[i]], [By[i], By[i] + vBy[i]])
    line_aB.set_data([Bx[i], Bx[i] + aBx[i]], [By[i], By[i] + aBy[i]])

    line_vC.set_data([Cx[i], Cx[i] + vCx[i]], [Cy[i], Cy[i] + vCy[i]])

    line_vD.set_data([Dx[i], Dx[i] + vDx[i]], [Dy[i], Dy[i] + vDy[i]])

    line_vE.set_data([Ex[i], Ex[i] + vEx[i]], [Ey[i], Ey[i] + vEy[i]])

    line_vF.set_data([Fx[i], Fx[i] + vFx[i]], [Fy[i], Fy[i] + vFy[i]])

    text.set_text("time = %.2f \n wAO1 = %.2f \n wBO2 = %.2f \n wAB = %.2f \n wAC = %.2f \n wBC = %.2f "
                  "\n wCE = %.2f \n wCD = %.2f \n wEF = %.2f \n wFO3 = %.2f"%
                  (t[i], wAO1[i], wBO2[i], wAB[i], wAC[i], wBC[i], wCE[i], wCD[i], wEF[i], wFO3[i]))
    return line, line_vA, line_aA, line_vB, line_aB, line_vC, line_vD, line_vE, line_vF, text

animation = FuncAnimation(fig, func=animation_frame, frames=len(t), interval=10, blit=True)
plt.vlines(np.arange(50, 110, 1), 15, 17, colors='black', linestyles='dashed')

# Uncomment to generate gif anumation
#savegif = PillowWriter(fps=30)
#animation.save("TM_HW2_Task1.gif", writer=savegif)

plt.show()