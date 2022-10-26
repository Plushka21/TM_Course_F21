import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Create plot
fig, ax = plt.subplots()
ax.set_xlim(-40, 80)
ax.set_ylim(-60, 60)
ax.set_aspect("equal")
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
lineM, = ax.plot([], [], 'o-', lw=2, label="M")
lineArc, = ax.plot([], [], lw=2)
line_vMtr, = ax.plot([], [], '<-', lw=2, label="v_M_tr")
line_vMrel, = ax.plot([], [], '<-', lw=2, label="v_M_rel")
line_vMabs, = ax.plot([], [], '<-', lw=2, label="v_M_abs")

line_aMtr, = ax.plot([], [], '<-', lw=2, label="a_M_tr")
line_aMrel, = ax.plot([], [], '<-', lw=2, label="a_M_rel")
line_aMabs, = ax.plot([], [], '<-', lw=2, label="a_M_abs")

# Text where angular velocity will be shown
text = ax.text(0.2, 0.2, 'matplotlib', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
leg = ax.legend(loc='upper left')

dt = 0.01
t = np.arange(0, np.pi, dt)
R = 18
k = 20

# Create required dataset
Ox, Oy = [], []
Ax, Ay = [], []
Mx, My = [], []
vMtr, vMxtr, vMytr, vMrel, vMxrel, vMyrel, vMabs, vMxabs, vMyabs = [], [], [], [], [], [], [], [], []
aMtr, aMxtr, aMytr, aMrel, aMxrel, aMyrel, aMabs, aMxabs, aMyabs = [], [], [], [], [], [], [], [], []
phi, al, arc = [], [], []
aO1, aO2 = [], []
bO1, bO2 = [], []

# Function to find positions of points O and A
def findOA():
    for i in range(len(t)):
        phi_cur = np.pi * t[i]**3 / 6
        phi.append(phi_cur)

        Ox.append(k * np.cos(phi_cur))
        Oy.append(k * np.sin(phi_cur))

        Ax.append(2 * R + k * np.cos(phi_cur))
        Ay.append(k * np.sin(phi_cur))

# Function to find position of point M
def findM():
    for i in range(len(t)):
        al_cur = (6 * np.pi * t[i] ** 2) / R

        if ((al_cur // np.pi) % 2 != 0):
            al_cur = -al_cur
        al.append(al_cur)

        Mx.append(Ox[i] + R * (1 - np.cos(al_cur)))
        My.append(Oy[i] + R * np.sin(al_cur))

# Function to find velicoty of point M
def findvM():
    for i in range(len(t)):
        om_cur = np.pi * t[i] ** 2 / 2 # omega(t) = phi(t)'
        vMtr_cur = k * om_cur
        vMtr.append(vMtr_cur)

        vMxtr.append(-vMtr_cur * np.cos(np.pi / 2 - phi[i]))
        vMytr.append( vMtr_cur * np.sin(np.pi / 2 - phi[i]))

        vMrel_cur = 12 * np.pi * t[i]  # v_rel(t) = s(t)'
        vMrel.append(vMrel_cur)

        vMxrel.append(vMrel_cur * np.cos(np.pi / 2 - al[i]))
        vMyrel.append(vMrel_cur * np.sin(np.pi / 2 - al[i]))
        if (al[i] < 0):
            vMxrel[i] = vMxrel[i] * (-1)
            vMyrel[i] = vMyrel[i] * (-1)

        vMxabs.append(vMxtr[i] + vMxrel[i])
        vMyabs.append(vMytr[i] + vMyrel[i])
        vMabs.append(np.sqrt(vMxabs[i] ** 2 + vMyabs[i] ** 2))

# Function to find acceleration of point M
def findaM():
    for i in range(len(t)):
        om_cur = np.pi * t[i] ** 2 / 2 # omega(t) = phi(t)'
        a_tr_n = k * om_cur ** 2
        eps_cur = np.pi * t[i] # eps(t) = omega(t)' = phi(t)''
        a_tr_t = k * eps_cur

        aMxtr.append(-a_tr_n * np.cos(phi[i]) - a_tr_t * np.cos(np.pi / 2 - phi[i]))
        aMytr.append(-a_tr_n * np.sin(phi[i]) + a_tr_t * np.sin(np.pi / 2 - phi[i]))
        aMtr.append(np.sqrt(aMxtr[i] ** 2 + aMytr[i] ** 2))

        a_rel_n = (12 * np.pi * t[i] ** 2) / R # (v_M)^2 / R
        a_rel_t = 12 * np.pi                   # (v_M)'

        if (al[i] > 0):
            aMxrel.append( a_rel_n * np.cos(al[i]) + a_rel_t * np.cos(np.pi / 2 - al[i]))
            aMyrel.append(-a_rel_n * np.sin(al[i]) + a_rel_t * np.sin(np.pi / 2 - al[i]))
        else:
            aMxrel.append( a_rel_n * np.cos(al[i]) - a_rel_t * np.cos(np.pi / 2 - al[i]))
            aMyrel.append(-a_rel_n * np.sin(al[i]) - a_rel_t * np.sin(np.pi / 2 - al[i]))

        aMrel.append(np.sqrt(aMxrel[i] ** 2 + aMyrel[i] ** 2))

        aMxabs.append(aMxtr[i] + aMxrel[i])
        aMyabs.append(aMytr[i] + aMyrel[i])
        aMabs.append(np.sqrt(aMxabs[i] ** 2 + aMyabs[i] ** 2))

findOA()
findM()
findvM()
findaM()
arc = [i for i in np.arange(0, np.pi, dt)]

# Generate animation
def animation_frame(i):
    line.set_data([0,     Ox[i], Ax[i], 36   ],
                  [Oy[0], Oy[i], Ay[i], Ay[0]])

    lineArc.set_data([R + Ox[i] + R * np.cos(arc)],
                     [Oy[i] + R * np.sin(arc)])

    lineM.set_data([Mx[i]], [My[i]])

    line_vMtr.set_data([Mx[i], Mx[i] + vMxtr[i]],
                       [My[i], My[i] + vMytr[i]])

    line_vMrel.set_data([Mx[i], Mx[i] + vMxrel[i]],
                        [My[i], My[i] + vMyrel[i]])

    line_vMabs.set_data([Mx[i], Mx[i] + vMxabs[i]],
                        [My[i], My[i] + vMyabs[i]])

    line_aMtr.set_data([Mx[i], Mx[i] + aMxtr[i]],
                       [My[i], My[i] + aMytr[i]])

    line_aMrel.set_data([Mx[i], Mx[i] + aMxrel[i]],
                        [My[i], My[i] + aMyrel[i]])

    line_aMabs.set_data([Mx[i], Mx[i] + aMxabs[i]],
                        [My[i], My[i] + aMyabs[i]])

    text.set_text("time = %.2f \n vM_tr = %.2f \n vM_rel = %.2f \n vM_abs = %.2f "
                  "\n aM_tr = %.2f \n aM_rel = %.2f \n aM_abs = %.2f"
                  %(t[i], vMtr[i], vMrel[i], vMabs[i], aMtr[i], aMrel[i], aMabs[i]))

    return line, lineM, lineArc, line_vMtr, line_vMrel, line_vMabs, line_aMtr, line_aMrel, line_aMabs, text

animation = FuncAnimation(fig, func=animation_frame, frames=len(t), interval=25, blit=True)

# Uncomment to generate gif anumation
#savegif = PillowWriter(fps=30)
#animation.save("TM_HW3_Task1.gif", writer=savegif)

plt.show()