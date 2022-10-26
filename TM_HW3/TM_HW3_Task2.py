import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Create plot
fig, ax = plt.subplots()
ax.set_xlim(-150, 100)
ax.set_ylim(-100, 100)
ax.set_aspect("equal")
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
lineArc, = ax.plot([], [], lw=2)
lineO1, = ax.plot([], [], 'o-', lw=2, label="O1")
lineO, = ax.plot([], [], 'o-', lw=2, label="O")
lineM, = ax.plot([], [], 'o-', lw=2, label="M")
line_vMtr, = ax.plot([], [], '<-', lw=2, label="v_M_tr")
line_vMrel, = ax.plot([], [], '<-', lw=2, label="v_M_rel")
line_vMabs, = ax.plot([], [], '<-', lw=2, label="v_M_abs")

line_aMtr, = ax.plot([], [], '<-', lw=2, label="a_M_tr")
line_aMrel, = ax.plot([], [], '<-', lw=2, label="a_M_rel")
line_aMabs, = ax.plot([], [], '<-', lw=2, label="a_M_abs")

# Text where angular velocity will be shown
text = ax.text(0.2, 0.2, 'matplotlib', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
leg = ax.legend(loc='upper left')

# Create required dataset
Cx, Cy = [], []
Ox, Oy = [], []
Mx, My = [], []
vM_tr, vMxtr, vMytr, vM_rel, vMxrel, vMyrel, vM_abs, vMxabs, vMyabs = [], [], [], [], [], [], [], [], []
aM_tr, aMxtr, aMytr, aM_rel, aMxrel, aMyrel, aM_abs, aMxabs, aMyabs = [], [], [], [], [], [], [], [], []
phi, al = [], []

dt = 0.01
t = np.arange(0, np.pi, dt)
R = 30

# Function to find positions of points O and C where C is center
def findCO():
    for i in range(len(t)):
        phi_cur = 2 * t[i] - 0.3 * t[i] ** 2
        phi.append(phi_cur)
        Cx.append(-R * np.sin(phi_cur))
        Cy.append(R * np.cos(phi_cur))

        Ox.append(Cx[i] + R * np.cos(phi[i]))
        Oy.append(Cy[i] + R * np.sin(phi[i]))

# Function to find position of point M
def findM():
    for i in range(len(t)):
        al_cur = (75 * np.pi * (0.1 * t[i] + 0.3 * t[i] ** 2)) / R  # Current angle OCm where C - center of circle
        k = al_cur // (2 * np.pi)                                   # Number of full turns of M
        al_cur = al_cur - (2 * k * np.pi)                           # After each full turn start counting alpha from 0
        al.append(al_cur)

        Mx.append(Cx[i] - R * np.cos(np.pi - al_cur - phi[i]))
        My.append(Cy[i] + R * np.sin(np.pi - al_cur - phi[i]))

# Function to find velocity of point M
def find_vM():
    temp_an = 3 * np.pi / 2
    for i in range(len(t)):
        vM_rel_cur = 75 * np.pi * (0.1 + 0.6 * t[i]) # = s(t)'
        vM_rel.append(vM_rel_cur)

        vMxrel.append(- vM_rel_cur * np.sin(np.pi - al[i] - phi[i]))
        vMyrel.append(- vM_rel_cur * np.cos(np.pi - al[i] - phi[i]))


        om_cur = 2 - 0.6 * t[i] # = phi(t)'
        MO1 = np.sqrt(2 * R ** 2 * (1 - np.cos(np.pi / 2 + al[i])))
        vM_tr_cur = om_cur * MO1
        vM_tr.append(vM_tr_cur)

        gamma1 = np.pi / 2 - (np.pi - (al[i] + phi[i]))
        gamma2 = np.pi / 4 - al[i] / 2
        gamma = gamma1 + gamma2
        vMxtr.append(- vM_tr_cur * np.cos(gamma))
        vMytr.append(- vM_tr_cur * np.sin(gamma))
        # When point M reaches point O1 direction of vector has to be inversed
        if ((al[i] // temp_an) % 2 != 0):
            vMxtr[i] = -vMxtr[i]
            vMytr[i] = -vMytr[i]

        vMxabs.append(vMxrel[i] + vMxtr[i])
        vMyabs.append(vMyrel[i] + vMytr[i])
        vM_abs.append(np.sqrt(vMxabs[i] ** 2 + vMyabs[i] ** 2))

# Function to find acceleration of point M
def find_aM():
    temp_an = 3 * np.pi / 2
    for i in range(len(t)):
        aM_rel_n = (75 * np.pi * (0.1 + 0.6 * t[i])) ** 2 / R   # = (vM)^2 / R
        aM_rel_t = 75 * np.pi * 0.6                             # = (vM)'

        beta = np.pi - al[i] - phi[i]
        aMxrel.append(aM_rel_n * np.cos(beta) - aM_rel_t * np.sin(beta))
        aMyrel.append(-aM_rel_n * np.sin(beta) - aM_rel_t * np.cos(beta))
        aM_rel.append(np.sqrt(aMxrel[i] ** 2 + aMyrel[i] ** 2))

        om_cur = 2 - 0.6 * t[i]    # = phi'
        MO1 = np.sqrt(2 * R ** 2 * (1 - np.cos(np.pi / 2 + al[i])))
        aM_tr_n = MO1 * om_cur ** 2
        eps_cur = -0.6             # = om_cur'
        aM_tr_t = MO1 * eps_cur

        gamma1 = np.pi / 2 - (np.pi - (al[i] + phi[i]))
        gamma2 = np.pi / 4 - al[i] / 2
        gamma = gamma1 + gamma2
        # When point M reaches point O1 direction of vector has to be inversed
        if ((al[i] // temp_an) % 2 != 0):
            aMxtr.append(- aM_tr_n * np.sin(gamma) + aM_tr_t * np.cos(gamma))
            aMytr.append(  aM_tr_n * np.cos(gamma) + aM_tr_t * np.sin(gamma))
        else:
            aMxtr.append(  aM_tr_n * np.sin(gamma) - aM_tr_t * np.cos(gamma))
            aMytr.append(- aM_tr_n * np.cos(gamma) - aM_tr_t * np.sin(gamma))

        aM_tr.append(np.sqrt(aMxtr[i] ** 2 + aMytr[i] ** 2))

        aMxabs.append(aMxtr[i] + aMxrel[i])
        aMyabs.append(aMytr[i] + aMyrel[i])
        aM_abs.append(np.sqrt(aMxabs[i] ** 2 + aMyabs[i] ** 2))


findCO()
findM()
find_vM()
find_aM()

arc = [i for i in np.arange(0, 2 * np.pi, dt)]

# Generate animation
def animation_frame(i):
    lineArc.set_data([Cx[i] + R * np.cos(arc)],
                     [Cy[i] + R * np.sin(arc)])

    lineO1.set_data([0], [0])

    lineO.set_data([Ox[i]],[Oy[i]])

    lineM.set_data([Mx[i]],[My[i]])

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
                  %(t[i], vM_tr[i], vM_rel[i], vM_abs[i], aM_tr[i], aM_rel[i], aM_abs[i]))
    return line, lineO, lineM, lineArc, line_vMtr, line_vMrel, line_vMabs, line_aMtr, line_aMrel, line_aMabs, text

animation = FuncAnimation(fig, func=animation_frame, frames=len(t), interval=25, blit=True)

# Uncomment to generate gif anumation
#savegif = PillowWriter(fps=30)
#animation.save("TM_HW3_Task2.gif", writer=savegif)

plt.show()