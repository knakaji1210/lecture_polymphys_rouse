# Rouse模型構築のためのステップ
# damped harmonic oscillatorからmassを外したバージョン

# ordinary differential equation of damped harmonic oscillator

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dampedHarmonicOscillator_nomass(x, t, c, k):

    v = -(k/c)*x

    return v

# variables
try:
    k = float(input('spring constant [N/m] (default=1.0): '))
except ValueError:
    k = 1.0                # [N/m] spring constant
try:
    c = float(input('damping coefficient [kg/s] (default=0.5): '))    # [kg/s] damping coefficient
except ValueError:
    c = 0.5                 # [m/s] damping coefficient

tau = c/k                   # [s] relaxation time

l = 20                      # [m] equilibrium length

tmax = 10*tau                # [s] duration time
dt = 0.05                   # [s] interval time

# initial condition
try:
    x0 = float(input('initial position [m] (default=10.0): '))
except ValueError:
    x0 = 10.0

t = np.arange(0, tmax, dt)

sol = odeint(dampedHarmonicOscillator_nomass, x0, t, args=(c,k))
x = sol[:, 0]

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-5, 50), ylim=(-5, 5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
var_template = r'$k$ = {0:.1f} N/m, $c$ = {1:.1f} kg/s'.format(k,c)
ax.text(0.4, 0.9, var_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
rtime_template = r'$\tau$ = {0:.2f} s'.format(tau)
ax.text(0.1, 0.8, rtime_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

line, = plt.plot([], [], 'ro-', animated=True)
# ここでは[],[]としているが、下でline.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():                 # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return line, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_pos = [0, l + x[i]]
    y_pos = [0, 0]
    line.set_data(x_pos,y_pos)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/damped_harmonicOsci_nomass_(x={0:.1f},k={1:.1f},c={2:.1f}).gif'.format(x0,k,c)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()