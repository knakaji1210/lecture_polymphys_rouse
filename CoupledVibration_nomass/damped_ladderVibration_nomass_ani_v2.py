# Rouse模型構築のためのステップ
# damped ladder vibrationからmassを外したバージョン

# ordinary differential equation of damped ladder vibration (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dampedLadderVibration_nomass(s, t, k, c):
    x1, x2, x3, x4, x5, x6, x7, x8, x9  = s

    v1 =             - (2*k/c) *x1 + (k/c) *x2
    v2 =   (k/c) *x1 - (2*k/c) *x2 + (k/c) *x3
    v3 =   (k/c) *x2 - (2*k/c) *x3 + (k/c) *x4
    v4 =   (k/c) *x3 - (2*k/c) *x4 + (k/c) *x5
    v5 =   (k/c) *x4 - (2*k/c) *x5 + (k/c) *x6
    v6 =   (k/c) *x5 - (2*k/c) *x6 + (k/c) *x7
    v7 =   (k/c) *x6 - (2*k/c) *x7 + (k/c) *x8
    v8 =   (k/c) *x7 - (2*k/c) *x8 + (k/c) *x9
    v9 =   (k/c) *x8 - (2*k/c) *x9

    dsdt = [v1, v2, v3, v4, v5, v6, v7, v8, v9]

    return dsdt

'''
上のladderVibrationをもっと簡潔に書けないかは未決
'''

# variables
try:
    k = float(input('spring constant [N/m] (default=1.0): '))
except ValueError:
    k = 1.0                # [N/m] spring constant
try:
    c = float(input('damping coefficient [kg/s] (default=0.5): '))    # [kg/s] damping coefficient
except ValueError:
    c = 0.5                 # [m/s] damping coefficient

l = 10                      # [m] equilibrium length for each component
n = 9                       # number of point
pos = [i*l for i in range(n+2)]
L = pos[-1]                 # [m] total length
tau = [(c/(4*k))/(np.sin((i+1)*np.pi/(2*(n+1))))**2 for i in range(n)]    # angular frequencies

try:
    mode2show = int(input('mode to visualize (default=4): '))
except ValueError:
    mode2show = 4

tmax = 3*np.max(tau)  # [s] duration time   

dt = 0.05                   # [s] interval time

# initial condition
amp = 5
try:
    ini = str(input('initial condition: "s"ingle, "m"ultiple : '))
except ValueError:
    ini = "s"
if ini == "s":
    try:
        mode = int(input('mode to stimulate (default=1): '))
    except ValueError:
        mode = 1
    s0 = [amp*np.sin(mode*np.pi*(i+1)/(n+1)) for i in range(n)]     # [m]   initial position
elif ini == "m":
    s0 = [amp*np.sin(np.pi*(i+1)/(n+1)) + amp*np.sin(3*np.pi*(i+1)/(n+1)) + amp*np.sin(9*np.pi*(i+1)/(n+1)) for i in range(n)]     # [m]   initial position
    mode = 9
else:
    s0 = np.zeros(n)     # [m]   initial position
    mode = 1

t = np.arange(0, tmax, dt)

sol = odeint(dampedLadderVibration_nomass, s0, t, args=(k,c))
x1, x2, x3, x4, x5, x6, x7, x8, x9 = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], sol[:, 4], sol[:, 5], sol[:, 6], sol[:, 7], sol[:, 8]
q1 = np.sin(np.pi*1/10)*x1 + np.sin(np.pi*2/10)*x2 + np.sin(np.pi*3/10)*x3 + np.sin(np.pi*4/10)*x4 + np.sin(np.pi*5/10)*x5 + np.sin(np.pi*6/10)*x6 + np.sin(np.pi*7/10)*x7 + np.sin(np.pi*8/10)*x8 + np.sin(np.pi*9/10)*x9                      # normal mode q1
q2 = np.sin(2*np.pi*1/10)*x1 + np.sin(2*np.pi*2/10)*x2 + np.sin(2*np.pi*3/10)*x3 + np.sin(2*np.pi*4/10)*x4 + np.sin(2*np.pi*5/10)*x5 + np.sin(2*np.pi*6/10)*x6 + np.sin(2*np.pi*7/10)*x7 + np.sin(2*np.pi*8/10)*x8 + np.sin(2*np.pi*9/10)*x9    # normal mode q2
q3 = np.sin(3*np.pi*1/10)*x1 + np.sin(3*np.pi*2/10)*x2 + np.sin(3*np.pi*3/10)*x3 + np.sin(3*np.pi*4/10)*x4 + np.sin(3*np.pi*5/10)*x5 + np.sin(3*np.pi*6/10)*x6 + np.sin(3*np.pi*7/10)*x7 + np.sin(3*np.pi*8/10)*x8 + np.sin(3*np.pi*9/10)*x9    # normal mode q3
qn = np.sin(mode2show*np.pi*1/10)*x1 + np.sin(mode2show*np.pi*2/10)*x2 + np.sin(mode2show*np.pi*3/10)*x3 + np.sin(mode2show*np.pi*4/10)*x4 + np.sin(mode2show*np.pi*5/10)*x5 + np.sin(mode2show*np.pi*6/10)*x6 + np.sin(mode2show*np.pi*7/10)*x7 + np.sin(mode2show*np.pi*8/10)*x8 + np.sin(mode2show*np.pi*9/10)*x9    # normal mode qn

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/5, 0.7*L))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
var_template = r'$k$ = {0:.1f} N/m, $c$ = {1:.1f} kg/s'.format(k,c)
ax.text(0.4, 0.92, var_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
rtime1_template = r'$\tau_1$ = {0:.2f} s, $\tau_2$ = {1:.2f} s, $\tau_3$ = {2:.2f} s'.format(tau[0],tau[1],tau[2])
ax.text(0.1, 0.85, rtime1_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
rtime2_template = r'$\tau_{0}$ = {1:.2f} s'.format(mode,tau[mode-1])
ax.text(0.1, 0.78, rtime2_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
norm3, = plt.plot([], [], 'yo-', animated=True)
norm, = plt.plot([], [], 'co-', animated=True)
# ここでは[],[]としているが、下でline.set_dataなどで実際の値を入れている

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.92, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return line, norm1, norm2, norm3, norm, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data(pos, 
            [0, x1[i], x2[i], x3[i], x4[i], x5[i], x6[i], x7[i], x8[i], x9[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [10, 10])
    norm2.set_data([0, L/2 + q2[i]], [20, 20])
    norm3.set_data([0, L/2 + q3[i]], [30, 30])
    norm.set_data([0, L/2 + qn[i]], [40, 40])
    time_text.set_text(time_template % (i*dt))
    return line, norm1, norm2, norm3, norm, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)

if ini == "s": 
    savefile = './gif/damped_ladderdVibration_single_(ini{0},show{1},k={2:.1f},c={3:.1f}).gif'.format(mode,mode2show,k,c)
    ani.save(savefile, writer='pillow', fps=fps)
elif ini == "m":
    savefile = './gif/damped_ladderdVibration_multiple_(show{0},k={1:.1f},c={2:.1f}).gif'.format(mode2show,k,c)
    ani.save(savefile, writer='pillow', fps=fps)
elif ini == "d":
    savefile = './gif/damped_ladderdVibration_delta_(show{0},k={1:.1f},c={2:.1f}).gif'.format(mode2show,k,c)
    ani.save(savefile, writer='pillow', fps=fps)
else:
    pass

plt.show()