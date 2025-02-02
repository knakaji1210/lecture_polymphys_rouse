# Rouse模型構築のためのステップ
# damped tripled harmonic oscillatorからmassを外したバージョン

# ordinary differential equation of damped tripled harmonic oscillator (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dampedTripledHarmonicOscillator_nomass(s, t, k, c):

    x1, x2, x3 = s              # s = (x1, x2, x3)
    v1 =             - (2*k/c) *x1 + (k/c) *x2
    v2 =   (k/c) *x1 - (2*k/c) *x2 + (k/c) *x3
    v3 =   (k/c) *x2 - (2*k/c) *x3
    dsdt = [v1, v2, v3]

    return dsdt

# variables
try:
    k = float(input('spring constant [N/m] (default=1.0): '))
except ValueError:
    k = 1.0               # [N/m] spring constant 2
try:
    c = float(input('damping coefficient [kg/s] (default=0.5): '))    # [kg/s] damping coefficient
except ValueError:
    c = 0.5                 # [m/s] damping coefficient

l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
l4 = 20                     # [m] equilibrium length
L = l1 + l2 + l3 + l4

tau1 = c/((2-np.sqrt(2))*k)     # [s] relaxation time
tau2 = c/(2*k)                  # [s] relaxation time
tau3 = c/((2+np.sqrt(2))*k)     # [s] relaxation time

tmax = 6*np.max([tau1,tau2,tau3])    # [s] duration time
dt = 0.05                       # [s] interval time


# initial condition
try:
    x1_0 = float(input('initial position of point1 (default=5.0): '))
except ValueError:
    x1_0 = 5.0
try:
    x2_0 = float(input('initial position of point2 (default=10.0): '))
except ValueError:
    x2_0 = 10.0
try:
    x3_0 = float(input('initial position of point3 (default=-5.0): '))
except ValueError:
    x3_0 = -5.0


s0 = [x1_0, x2_0, x3_0]   # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(dampedTripledHarmonicOscillator_nomass, s0, t, args=(k, c))  # ODEの解を求めている
x1, x2, x3 = sol[:, 0], sol[:, 1], sol[:, 2]    # [x1], [x2], [x3]が出てくる
q1 = (1/2)*(x1 + np.sqrt(2)*x2+x3)              # normal mode q1
q2 = (1/2)*(-np.sqrt(2)*x1 + np.sqrt(2)*x3)     # normal mode q2
q3 = (1/2)*(x1 - np.sqrt(2)*x2+x3)              # normal mode q3

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/4, 1.4*L))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
var_template = r'$k$ = {0:.1f} N/m, $c$ = {1:.1f} kg/s'.format(k,c)
ax.text(0.4, 0.92, var_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
rtime_template = r'$\tau_1$ = {0:.2f} s, $\tau_2$ = {1:.2f} s, $\tau_3$ = {2:.2f} s'.format(tau1,tau2,tau3)
ax.text(0.1, 0.78, rtime_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
norm3, = plt.plot([], [], 'yo-', animated=True)
# ここでは[],[]としているが、下でlinei.set_dataで実際の値を入れている

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.92, '', transform=ax.transAxes)
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return line, norm1, norm2, norm3, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l1, l1 + l2, l1 + l2 + l3, L], [0, x1[i], x2[i], x3[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [l1, l1])
    norm2.set_data([0, L/2 + q2[i]], [2*l1, 2*l1])
    norm3.set_data([0, L/2 + q3[i]], [3*l1, 3*l1])
    time_text.set_text(time_template % (i*dt))
    return line, norm1, norm2, norm3, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/damped_tripledHarmonicOsci_t_(x1={0:.1f},x2={1:.1f},x3={2:.1f},k={3:.1f},c={4:.1f}).gif'.format(x1_0,x2_0,x3_0,k,c)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()