#!/usr/bin/env python3

import math

import sympy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as pwd

from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider

spd = 40.0

t = sympy.symbols('t')

def bezier_sym(n) :
	p_lst = sympy.symbols(' '.join(f'P{i}' for i in range(n)))

	c_lst = [1, 1]
	for i in range(n-2) :
		c_lst = [1,] + [a + b for a, b in zip(c_lst[1:], c_lst[:-1])] + [1,]

	res = 0
	for i, (c, p) in enumerate(zip(c_lst, p_lst)) :
		res += t**(i) * (1-t)**(n-i-1) * c * p
	print(c_lst, p_lst)

	return res.expand().collect('t')

def b9func(t, * P) :
	assert(len(P) == 9)
	return P[0] + t**8*(P[0] - 8*P[1] + 28*P[2] - 56*P[3] + 70*P[4] - 56*P[5] + 28*P[6] - 8*P[7] + P[8]) + t**7*(-8*P[0] + 56*P[1] - 168*P[2] + 280*P[3] - 280*P[4] + 168*P[5] - 56*P[6] + 8*P[7]) + t**6*(28*P[0] - 168*P[1] + 420*P[2] - 560*P[3] + 420*P[4] - 168*P[5] + 28*P[6]) + t**5*(-56*P[0] + 280*P[1] - 560*P[2] + 560*P[3] - 280*P[4] + 56*P[5]) + t**4*(70*P[0] - 280*P[1] + 420*P[2] - 280*P[3] + 70*P[4]) + t**3*(-56*P[0] + 168*P[1] - 168*P[2] + 56*P[3]) + t**2*(28*P[0] - 56*P[1] + 28*P[2]) + t*(-8*P[0] + 8*P[1])

def b10int(t, * P) :
	assert(len(P) == 9)
	return P[0]*t + t**9*(P[0]/9 - 8*P[1]/9 + 28*P[2]/9 - 56*P[3]/9 + 70*P[4]/9 - 56*P[5]/9 + 28*P[6]/9 - 8*P[7]/9 + P[8]/9) + t**8*(-P[0] + 7*P[1] - 21*P[2] + 35*P[3] - 35*P[4] + 21*P[5] - 7*P[6] + P[7]) + t**7*(4*P[0] - 24*P[1] + 60*P[2] - 80*P[3] + 60*P[4] - 24*P[5] + 4*P[6]) + t**6*(-28*P[0]/3 + 140*P[1]/3 - 280*P[2]/3 + 280*P[3]/3 - 140*P[4]/3 + 28*P[5]/3) + t**5*(14*P[0] - 56*P[1] + 84*P[2] - 56*P[3] + 14*P[4]) + t**4*(-14*P[0] + 42*P[1] - 42*P[2] + 14*P[3]) + t**3*(28*P[0]/3 - 56*P[1]/3 + 28*P[2]/3) + t**2*(-4*P[0] + 4*P[1])

p_lst = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]

class BezierExp() :
	def __init__(self, speed, coefw=0.1) :
		self.speed = speed
		self.coefw = coefw

		self.ct = CoordinatedTurn(12.0, 30.0)
		self.dg = DummyGlider(lon=420.0, vx=spd)

		self.p_lst = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]

		self.line_map = dict()

	def setup(self) :
		self.fig = plt.figure()

		plt.subplot(3,3,1)
		self.line_map['control'], = plt.plot([90.0 * i/(8 * b10_final) for i in range(9)], np.array(p_lst), 'o-', linewidth=5)
		self.line_map['spline'], = plt.plot(u_arr, b9func(t_arr, * p_lst))
		plt.grid()


	def update(self) :



b10_final = b10int(1.0, * p_lst)
k = 90.0 / b10_final
j = 0.1

t_arr = np.linspace(0.0, 1.0, 1000)

for i in range(80000) :
	d = j * dg.lon / dg.vx
	psi = -k * b10int(min(1.0, abs(d)), * p_lst)
	if i == 0 :
		print(f"d={d} psi={psi}")
	dg.step_vx_psi(psi=psi)
	if dg.lon <= 0.25 and abs(dg.psi) <= 1.0 :
		break

m = dg.freeze()

u_arr = k * t_arr

plt.subplot(3,3,1)
line_control_points, = plt.plot([90.0 * i/(8 * b10_final) for i in range(9)], np.array(p_lst), 'o-', linewidth=5)
plt.plot(u_arr, b9func(t_arr, * p_lst))
plt.grid()

print(line_control_points)

plt.subplot(3,3,4)
plt.plot(u_arr, k * b10int(t_arr, * p_lst))
plt.grid()

plt.subplot(1, 3, 2)
plt.plot(m['lon'], m['lat'])
plt.title("lat / lon")
plt.grid()
plt.axis('equal')

plt.subplot(3, 3, 3)
plt.plot(m['t'], m['psi'])
plt.grid()
plt.subplot(3, 3, 6)
m['psidot'] = np.hstack(([math.nan,], (m['psi'][1:] - m['psi'][:-1]) / dg.dt))
plt.plot(m['t'], m['psidot'])
plt.axhline(ct.psidot_deg(3600 * spd / 1852), color='tab:red')
plt.grid()
plt.subplot(3, 3, 9)
m['psiddot'] = np.hstack(([math.nan,], (m['psidot'][1:] - m['psidot'][:-1]) / dg.dt))
plt.plot(m['t'], m['psiddot'])
plt.grid()

a_lst = [
	fig.add_axes([0.05 + 0.032*i, 0.1, 0.02, 0.20]) for i in range(9)
]
s_lst = [
	pwd.Slider(ax=a_lst[i], label=f"P{i}", valmin=-0.333, valmax=1.333, valinit=p_lst[i], orientation="vertical") for i in range(9)
]


def update(val) :
	p_lst = [s.val for s in s_lst]

	line_control_points.set_ydata(p_lst)
	fig.canvas.draw_idle()

for s in s_lst :
	s.on_changed(update)

plt.show()

# if __name__ == '__main__' :
# 	print(bezier_sym(4))
# 	b9 = bezier_sym(9)
# 	print(b9)
# 	b10 = sympy.integrate(b9, t)
# 	print(b10)
	
