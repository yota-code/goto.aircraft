#!/usr/bin/env python3

# https://fr.wikipedia.org/wiki/Courbe_du_chien
# https://www.tutorialspoint.com/signals-and-systems-what-is-inverse-z-transform#
# https://lpsa.swarthmore.edu/ZXform/InvZXform/InvZXform.html

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

# if True :
# 	b10fnc = bezier_sym(10)
# 	print(b10fnc)
# 	b11int = sympy.integrate(b10fnc, t)
# 	print(b11int)
# 	sys.exit(0)

p_lst = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]

class BezierPlot() :

	def __init__(self, n, speed=30.0, coefw=0.1) :
		self.n = n

		self.speed = speed
		self.coefw = coefw
		self.gamma = 1.1

		self.ct = CoordinatedTurn(12.0, 30.0)

		if n == 9 :
			self.p_lst = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]
			self.p_lst = [1.0, 1.0, 0.094, 0.745, 0.215, 0.559, 0.301, 0.0, 0.0]
			self.p_lst = [1.0, 1.0, 1.0, 0.044, 0.616, 1.0, 0.0, 0.0, 0.0]
		elif n == 10 :
			# self.p_lst = [1.0, 1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]
			self.p_lst = [1.0, 1.0, 1.0, 0.05, 0.916, 0.137, 0.96, 0.0, 0.0, 0.0]
			# self.p_lst = [1.0, 1.0, 1.0, 0.05, 0.916, 0.137, 0.273, 0.087, 0.0, 0.0]

		self.line_map = dict()
		self.axe_map = dict()
		self.score_map = dict()

		self.plot_setup()
		self.plot_update(0.0)

	def score(self, p_lst) :
		p_tup = tuple(p_lst)

		psidot_max = self.ct.psidot_deg(3600 * spd / 1852)

		if p_tup in self.score_map :
			return self.score_map[p_tup]

		m = self.run(p_lst)

	def update(self, p_lst) :

		assert(len(p_lst) == self.n)

		self.p_lst = p_lst

		if self.n == 9 :
			self.c_lst = [
				self.p_lst[0], # [0,]
				-8*self.p_lst[0] + 8*self.p_lst[1],
				28*self.p_lst[0] - 56*self.p_lst[1] + 28*self.p_lst[2],
				-56*self.p_lst[0] + 168*self.p_lst[1] - 168*self.p_lst[2] + 56*self.p_lst[3],
				70*self.p_lst[0] - 280*self.p_lst[1] + 420*self.p_lst[2] - 280*self.p_lst[3] + 70*self.p_lst[4],
				-56*self.p_lst[0] + 280*self.p_lst[1] - 560*self.p_lst[2] + 560*self.p_lst[3] - 280*self.p_lst[4] + 56*self.p_lst[5],
				28*self.p_lst[0] - 168*self.p_lst[1] + 420*self.p_lst[2] - 560*self.p_lst[3] + 420*self.p_lst[4] - 168*self.p_lst[5] + 28*self.p_lst[6],
				-8*self.p_lst[0] + 56*self.p_lst[1] - 168*self.p_lst[2] + 280*self.p_lst[3] - 280*self.p_lst[4] + 168*self.p_lst[5] - 56*self.p_lst[6] + 8*self.p_lst[7],
				self.p_lst[0] - 8*self.p_lst[1] + 28*self.p_lst[2] - 56*self.p_lst[3] + 70*self.p_lst[4] - 56*self.p_lst[5] + 28*self.p_lst[6] - 8*self.p_lst[7] + self.p_lst[8]
			]
		elif self.n == 10 :
			self.c_lst = [
				self.p_lst[0],
				-9*self.p_lst[0] + 9*self.p_lst[1],
				36*self.p_lst[0] - 72*self.p_lst[1] + 36*self.p_lst[2],
				-84*self.p_lst[0] + 252*self.p_lst[1] - 252*self.p_lst[2] + 84*self.p_lst[3],
				126*self.p_lst[0] - 504*self.p_lst[1] + 756*self.p_lst[2] - 504*self.p_lst[3] + 126*self.p_lst[4],
				-126*self.p_lst[0] + 630*self.p_lst[1] - 1260*self.p_lst[2] + 1260*self.p_lst[3] - 630*self.p_lst[4] + 126*self.p_lst[5],
				84*self.p_lst[0] - 504*self.p_lst[1] + 1260*self.p_lst[2] - 1680*self.p_lst[3] + 1260*self.p_lst[4] - 504*self.p_lst[5] + 84*self.p_lst[6],
				-36*self.p_lst[0] + 252*self.p_lst[1] - 756*self.p_lst[2] + 1260*self.p_lst[3] - 1260*self.p_lst[4] + 756*self.p_lst[5] - 252*self.p_lst[6] + 36*self.p_lst[7],
				9*self.p_lst[0] - 72*self.p_lst[1] + 252*self.p_lst[2] - 504*self.p_lst[3] + 630*self.p_lst[4] - 504*self.p_lst[5] + 252*self.p_lst[6] - 72*self.p_lst[7] + 9*self.p_lst[8],
				-self.p_lst[0] + 9*self.p_lst[1] - 36*self.p_lst[2] + 84*self.p_lst[3] - 126*self.p_lst[4] + 126*self.p_lst[5] - 84*self.p_lst[6] + 36*self.p_lst[7] - 9*self.p_lst[8] + self.p_lst[9],
			]

		print(self.c_lst)

		self.final = self.b_inte(1.0)

	def run(self, p_lst) :

		self.update(p_lst)

		k = 90.0 / self.final

		self.dg = DummyGlider(lon=420.0, vx=self.speed)

		circle_radius = self.ct.circle_radius(self.speed)
		for i in range(80000) :
			d = max(-1.0, min(self.dg.lon / circle_radius, 1.0))
			psi = -k * self.b_inte(min(1.0, abs(d)))
			if i == 0 :
				print(f"d={d} psi={psi} p_lst={p_lst}")
			self.dg.step_vx_psi(psi=psi)
			if self.dg.lon <= 0.25 and abs(self.dg.psi) <= 1.0 :
				break

		return self.dg.freeze()

	def plot_setup(self) :
		self.fig = plt.figure()

		self.axe_map["3.4.1"] = plt.subplot(3, 3, 1)
		self.line_map['control'], = plt.plot([0.0,], 'o-', linewidth=5)
		self.line_map['b_func'], = plt.plot([0.0,])
		plt.grid()

		self.axe_map["3.4.5"] = plt.subplot(3, 3, 4)
		self.line_map['b_inte'], =plt.plot([0.0,])
		plt.grid()

		self.axe_map["1.4.2"] = plt.subplot(1, 3, 2)
		self.line_map['trace'], = plt.plot([0.0,])
		plt.title("lat / lon")
		plt.grid()
		plt.axis('equal')

		self.axe_map["3.4.3"] = plt.subplot(4, 3, 3)
		self.line_map['psi'], = plt.plot([0.0,])
		plt.ylabel('psi')
		plt.grid()
		self.axe_map["3.4.7"] = plt.subplot(4, 3, 6)
		self.line_map['psi1'], = plt.plot([0.0,])
		plt.axhline(self.ct.psidot_max_deg(v_ms=self.speed), color='tab:red')
		plt.ylabel('psidot (1)')
		plt.grid()
		self.axe_map["3.4.9"] = plt.subplot(4, 3, 9)
		self.line_map['psi2'], = plt.plot([0.0,])
		plt.axhline(self.ct.psiddot_max_deg(v_ms=self.speed), color='tab:red')
		plt.ylabel('psiddot (2)')
		plt.grid()
		self.axe_map["3.4.12"] = plt.subplot(4, 3, 12)
		self.line_map['psi3'], = plt.plot([0.0,])
		plt.ylabel('psidddot (3)')
		plt.grid()

		self.a_lst = [
			self.fig.add_axes([0.05 + 0.032*i, 0.1, 0.02, 0.20]) for i in range(self.n)
		]
		self.s_lst = [
			pwd.Slider(ax=self.a_lst[i], label=f"P{i}", valmin=-0.2, valmax=1.2, valinit=self.p_lst[i], orientation="vertical") for i in range(self.n)
		]
		for s in self.s_lst :
			s.on_changed(self.plot_update)

	def plot_update(self, val) :
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

		# self.final = self.b_inte(1.0)
		k = 90.0 / self.final
		u_arr = k * t_arr

		# self.dg = DummyGlider(lon=420.0, vx=self.speed)
		# for i in range(80000) :
		# 	d = self.coefw * self.dg.lon / self.dg.vx
		# 	psi = -k * self.b_inte(min(1.0, abs(d)))
		# 	if i == 0 :
		# 		print(f"d={d} psi={psi} p_lst={p_lst}")
		# 	self.dg.step_vx_psi(psi=psi)
		# 	if self.dg.lon <= 0.25 and abs(self.dg.psi) <= 1.0 :
		# 		break


		self.line_map['control'].set_xdata([90.0 * i/((self.n - 1) * self.final) for i in range(self.n)])
		self.line_map['control'].set_ydata(self.p_lst)

		self.line_map['b_func'].set_xdata(u_arr)
		self.line_map['b_func'].set_ydata(self.b_func(t_arr))

		self.line_map['b_inte'].set_xdata(u_arr)
		self.line_map['b_inte'].set_ydata(k * self.b_inte(t_arr))

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

		print("updated", self.p_lst)

		self.fig.canvas.draw_idle()

	def adapt(self, t) :
		return max(0.0, min(t, 1.0))**self.gamma

	def b_func(self, t) :
		return sum([self.c_lst[i]*self.adapt(t)**i for i in range(self.n)])

	def b_inte(self, t) :
		return sum([self.c_lst[i]*self.adapt(t)**(i+1)/(i+1) for i in range(self.n)])

u = BezierPlot(10)
plt.show()
