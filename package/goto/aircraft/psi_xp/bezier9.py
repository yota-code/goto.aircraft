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

p_lst = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]

class Bezier9plot() :
	def __init__(self, speed=30.0, coefw=0.1) :
		self.speed = speed
		self.coefw = coefw

		self.ct = CoordinatedTurn(12.0, 30.0)

		self.p_lst = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]
		self.p_lst = [1.0, 1.0, 0.094, 0.745, 0.215, 0.559, 0.301, 0.0, 0.0]
		self.p_lst = [1.0, 1.0, 1.0, 0.044, 0.616, 1.0, 0.0, 0.0, 0.0]

		self.line_map = dict()
		self.axe_map = dict()
		self.score_map = dict()

		self.setup()
		self.update(0.0)

	def score(self, p_lst) :
		p_tup = tuple(p_lst)

		psidot_max = self.ct.psidot_deg(3600 * spd / 1852)

		if p_tup in self.score_map :
			return self.score_map[p_tup]

		m = self.run(p_lst)


	def run(self, p_lst) :

		assert(len(p_lst) == 9)

		self.p_lst = p_lst
		self.c_lst = [
			self.p_lst[0],
			-8*self.p_lst[0] + 8*self.p_lst[1],
			28*self.p_lst[0] - 56*self.p_lst[1] + 28*self.p_lst[2],
			-56*self.p_lst[0] + 168*self.p_lst[1] - 168*self.p_lst[2] + 56*self.p_lst[3],
			70*self.p_lst[0] - 280*self.p_lst[1] + 420*self.p_lst[2] - 280*self.p_lst[3] + 70*self.p_lst[4],
			-56*self.p_lst[0] + 280*self.p_lst[1] - 560*self.p_lst[2] + 560*self.p_lst[3] - 280*self.p_lst[4] + 56*self.p_lst[5],
			28*self.p_lst[0] - 168*self.p_lst[1] + 420*self.p_lst[2] - 560*self.p_lst[3] + 420*self.p_lst[4] - 168*self.p_lst[5] + 28*self.p_lst[6],
			-8*self.p_lst[0] + 56*self.p_lst[1] - 168*self.p_lst[2] + 280*self.p_lst[3] - 280*self.p_lst[4] + 168*self.p_lst[5] - 56*self.p_lst[6] + 8*self.p_lst[7],
			self.p_lst[0] - 8*self.p_lst[1] + 28*self.p_lst[2] - 56*self.p_lst[3] + 70*self.p_lst[4] - 56*self.p_lst[5] + 28*self.p_lst[6] - 8*self.p_lst[7] + self.p_lst[8]
		]

		self.final = self.b10int(1.0)
		k = 90.0 / self.final

		print(self.p_lst)
		print(self.final, k)

		self.dg = DummyGlider(lon=420.0, vx=self.speed)

		for i in range(80000) :
			d = self.coefw * self.dg.lon / self.dg.vx
			psi = -k * self.b10int(min(1.0, abs(d)))
			if i == 0 :
				print(f"d={d} psi={psi} p_lst={p_lst}")
			self.dg.step_vx_psi(psi=psi)
			if self.dg.lon <= 0.25 and abs(self.dg.psi) <= 1.0 :
				break

		return self.dg.freeze()

	def setup(self) :
		self.fig = plt.figure()

		self.axe_map["3.4.1"] = plt.subplot(3, 3, 1)
		self.line_map['control'], = plt.plot([0.0,], 'o-', linewidth=5)
		self.line_map['b9func'], = plt.plot([0.0,])
		plt.grid()

		self.axe_map["3.4.5"] = plt.subplot(3, 3, 4)
		self.line_map['b10int'], =plt.plot([0.0,])
		plt.grid()

		self.axe_map["1.4.2"] = plt.subplot(1, 3, 2)
		self.line_map['trace'], = plt.plot([0.0,])
		plt.title("lat / lon")
		plt.grid()
		plt.axis('equal')

		self.axe_map["3.4.3"] = plt.subplot(4, 3, 3)
		self.line_map['psi'], = plt.plot([0.0,])
		plt.grid()
		self.axe_map["3.4.7"] = plt.subplot(4, 3, 6)
		self.line_map['psi1'], = plt.plot([0.0,])
		plt.axhline(self.ct.psidot_max_deg(v_ms=self.speed), color='tab:red')
		plt.grid()
		self.axe_map["3.4.9"] = plt.subplot(4, 3, 9)
		self.line_map['psi2'], = plt.plot([0.0,])
		plt.axhline(self.ct.psiddot_max_deg(v_ms=self.speed), color='tab:red')
		plt.grid()
		self.axe_map["3.4.12"] = plt.subplot(4, 3, 12)
		self.line_map['psi3'], = plt.plot([0.0,])
		plt.grid()

		self.a_lst = [
			self.fig.add_axes([0.05 + 0.032*i, 0.1, 0.02, 0.20]) for i in range(9)
		]
		self.s_lst = [
			pwd.Slider(ax=self.a_lst[i], label=f"P{i}", valmin=-0.2, valmax=1.2, valinit=self.p_lst[i], orientation="vertical") for i in range(9)
		]
		for s in self.s_lst :
			s.on_changed(self.update)

	def update(self, val) :
		t_arr = np.linspace(0.0, 1.0, 1000)

		m = self.run([s.val for s in self.s_lst])

		# self.p_lst = [s.val for s in self.s_lst]
		# self.c_lst = [
		# 	self.p_lst[0],
		# 	-8*self.p_lst[0] + 8*self.p_lst[1],
		# 	28*self.p_lst[0] - 56*self.p_lst[1] + 28*self.p_lst[2],
		# 	-56*self.p_lst[0] + 168*self.p_lst[1] - 168*self.p_lst[2] + 56*self.p_lst[3],
		# 	70*self.p_lst[0] - 280*self.p_lst[1] + 420*self.p_lst[2] - 280*self.p_lst[3] + 70*self.p_lst[4],
		# 	-56*self.p_lst[0] + 280*self.p_lst[1] - 560*self.p_lst[2] + 560*self.p_lst[3] - 280*self.p_lst[4] + 56*self.p_lst[5],
		# 	28*self.p_lst[0] - 168*self.p_lst[1] + 420*self.p_lst[2] - 560*self.p_lst[3] + 420*self.p_lst[4] - 168*self.p_lst[5] + 28*self.p_lst[6],
		# 	-8*self.p_lst[0] + 56*self.p_lst[1] - 168*self.p_lst[2] + 280*self.p_lst[3] - 280*self.p_lst[4] + 168*self.p_lst[5] - 56*self.p_lst[6] + 8*self.p_lst[7],
		# 	self.p_lst[0] - 8*self.p_lst[1] + 28*self.p_lst[2] - 56*self.p_lst[3] + 70*self.p_lst[4] - 56*self.p_lst[5] + 28*self.p_lst[6] - 8*self.p_lst[7] + self.p_lst[8]
		# ]

		# self.final = self.b10int(1.0)
		k = 90.0 / self.final
		u_arr = k * t_arr

		# self.dg = DummyGlider(lon=420.0, vx=self.speed)
		# for i in range(80000) :
		# 	d = self.coefw * self.dg.lon / self.dg.vx
		# 	psi = -k * self.b10int(min(1.0, abs(d)))
		# 	if i == 0 :
		# 		print(f"d={d} psi={psi} p_lst={p_lst}")
		# 	self.dg.step_vx_psi(psi=psi)
		# 	if self.dg.lon <= 0.25 and abs(self.dg.psi) <= 1.0 :
		# 		break


		self.line_map['control'].set_xdata([90.0 * i/(8 * self.final) for i in range(9)])
		self.line_map['control'].set_ydata(self.p_lst)

		self.line_map['b9func'].set_xdata(u_arr)
		self.line_map['b9func'].set_ydata(self.b9func(t_arr))

		self.line_map['b10int'].set_xdata(u_arr)
		self.line_map['b10int'].set_ydata(k * self.b10int(t_arr))

		self.line_map['trace'].set_xdata(m['lon'])
		self.line_map['trace'].set_ydata(m['lat'])

		m['psi1'] = np.hstack(([math.nan,], (m['psi'][1:] - m['psi'][:-1]) / self.dg.dt))
		m['psi2'] = np.hstack(([math.nan,], (m['psi1'][1:] - m['psi1'][:-1]) / self.dg.dt))
		m['psi3'] = np.hstack(([math.nan,], (m['psi2'][1:] - m['psi2'][:-1]) / self.dg.dt))

		for c in ['psi', 'psi1', 'psi2', 'psi3'] :
			self.line_map[c].set_xdata(m['t'])
			self.line_map[c].set_ydata(m[c])

		for k, v in self.axe_map.items() :
			v.relim()
			v.autoscale_view()

		print("updated", p_lst)

		self.fig.canvas.draw_idle()

	def b9func(self, t) :
		assert(len(self.p_lst) == 9)
		return self.c_lst[8]*t**8 + self.c_lst[7]*t**7 + self.c_lst[6]*t**6 + self.c_lst[5]*t**5 + self.c_lst[4]*t**4 + self.c_lst[3]*t**3 + self.c_lst[2]*t**2 + self.c_lst[1]*t + self.c_lst[0]

	def b10int(self, t) :
		assert(len(self.p_lst) == 9)
		return self.c_lst[8]*t**9/9 + self.c_lst[7]*t**8/8 + self.c_lst[6]*t**7/7 + self.c_lst[5]*t**6/6 + self.c_lst[4]*t**5/5 + self.c_lst[3]*t**4/4 + self.c_lst[2]*t**3/3 + self.c_lst[1]*t**2/2 + self.c_lst[0]*t



u = Bezier9plot()
plt.show()
sys.exit(0)



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
	
