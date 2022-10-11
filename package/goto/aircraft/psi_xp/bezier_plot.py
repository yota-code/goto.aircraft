#!/usr/bin/env python3

# https://fr.wikipedia.org/wiki/Courbe_du_chien
# https://www.tutorialspoint.com/signals-and-systems-what-is-inverse-z-transform#
# https://lpsa.swarthmore.edu/ZXform/InvZXform/InvZXform.html

import math
import sys

import sympy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as pwd

from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider

""" gamma met la zouille dans la version intégrale et pour l'instant on ne veut pas se passer de l'intégration donc j'ai enlevé gamma """

spd = 40.0

t = sympy.symbols('t')

def bezier_sym(n) :
	p_lst = sympy.symbols(' '.join(f'P{i}' for i in range(n)))
	g = sympy.symbols('gamma')

	c_lst = [1, 1]
	for i in range(n-2) :
		c_lst = [1,] + [a + b for a, b in zip(c_lst[1:], c_lst[:-1])] + [1,]

	res = 0
	for i, (c, p) in enumerate(zip(c_lst, p_lst)) :
		res += (t**i) * (1-t)**(n-i-1) * c * p

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

	def __init__(self, n, speed=30.0) :
		self.n = n
		self.n_lst = [i / (self.n-1) for i in range(self.n)]

		self.speed = speed
		self.slope = 1.0

		self.ct = CoordinatedTurn(12.0, 30.0)
		self.ct.display_info(self.speed)

		if n == 9 :
			self.p_lst = [1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]
			self.p_lst = [1.0, 1.0, 0.094, 0.745, 0.215, 0.559, 0.301, 0.0, 0.0]
			self.p_lst = [1.0, 1.0, 1.0, 0.044, 0.616, 1.0, 0.0, 0.0, 0.0]
		elif n == 10 :
			# self.p_lst = self.n_lst[::-1]
			# self.p_lst = [1.0 for i in self.n_lst]
			# self.p_lst = [1.0, 1.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0]
			# self.p_lst = [1.0, 1.0, 1.0, 0.1, 0.9, 0.1, 0.9, 0.0, 0.0, 0.0]
			# self.p_lst = [1.0, 1.0, 1.0, 0.05, 0.916, 0.137, 0.96, 0.0, 0.0, 0.0]
			# self.p_lst = [1.0, 1.0, 1.0, 0.05, 0.916, 0.137, 0.273, 0.087, 0.0, 0.0]
			# self.p_lst = [1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0]
			self.p_lst = [1.0, 1.0, 1.0, 0.08947368421052626, 0.8526315789473682, 0.2526315789473683, 0.8999999999999997, 0.0, 0.0, 0.0]
			self.p_lst = [1, 1, 0.152, 0.245, 0.812, 0.13, 0.467, 0, 0, 0]
			self.p_lst = [1, 0.93, 0.0367, 0.28, 0.945, 0.0, 0.18, 0.0, 0, 0]
			self.p_lst = [1, 0.909, 0.0367, 0.416, 0.816, 0, 0.18, 0, 0, 0]

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

		self.p_lst = [self.slope * i for i in p_lst]

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

		self.final = self.b_inte(1.0)

	def run(self) :

		reference_radius = self.ct.circle_radius(30.0)
		circle_radius = self.ct.circle_radius(self.speed)

		j = 90.0 / self.final
		
		self.dg = DummyGlider(lon=circle_radius * j * 1.2 / reference_radius, vx=self.speed)

		for i in range(80000) :
			d = reference_radius * self.dg.lon / (j * circle_radius)
			psi = - j * self.b_inte(min(1.0, abs(d)))
			if i == 0 :
				print("-" * 16)
				print("\tp_lst : [" + ', '.join(f"{i:.3g}" for i in self.p_lst) + "]")
				print("\tc_lst : [" + ', '.join(f"{i:.3g}" for i in self.c_lst) + "]")
				print(f"\tradius : {circle_radius}")
				print(f"\tfinal : {self.final}")
				print(f"\tpsi : {psi:.4g} (d = {d:.4g})")
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
		self.line_map['b_inte'], = plt.plot([0.0,])
		self.line_map['slope'], = plt.plot([0.0,])
		plt.grid()

		self.axe_map["1.4.2"] = plt.subplot(2, 3, 2)
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
		self.line_map['psidot_max'], = plt.plot([0.0,], color='tab:red')
		plt.ylabel('psidot (1)')
		plt.grid()
		self.axe_map["3.4.9"] = plt.subplot(4, 3, 9)
		self.line_map['psi2'], = plt.plot([0.0,])
		self.line_map['psiddot_max'], = plt.plot([0.0,], color='tab:red')
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

		self.s_map = dict()

		a_tmp = self.fig.add_axes([0.40, 0.36, 0.22, 0.02])
		self.s_map['speed'] = pwd.Slider(ax=a_tmp, label=f"speed", valmin=0.0, valmax=75.0, valinit=self.speed, orientation="horizontal")
		self.s_map['speed'].on_changed(self.plot_update)

		a_tmp = self.fig.add_axes([0.40, 0.4, 0.22, 0.02])
		self.s_map['slope'] = pwd.Slider(ax=a_tmp, label=f"slope", valmin=0.0, valmax=2.0, valinit=self.slope, orientation="horizontal")
		self.s_map['slope'].on_changed(self.plot_update)

	def plot_update(self, val) :
		u_arr = np.linspace(0.0, 1.0, 512)

		s_lst = [s.val for s in self.s_lst]

		self.speed = self.s_map['speed'].val
		self.slope = self.s_map['slope'].val

		self.update(s_lst)

		m = self.run()

		# self.line_map['control'].set_xdata(90.0 * ( np.linspace(0.0, 1.0, self.n)**(1.0 / self.gamma) ) / self.final)
		self.line_map['control'].set_xdata(self.n_lst)
		self.line_map['control'].set_ydata(self.p_lst)

		self.line_map['b_func'].set_xdata(u_arr)
		self.line_map['b_func'].set_ydata(self.b_func(u_arr))

		j = 90.0 / self.final
		self.line_map['b_inte'].set_xdata(j * u_arr)
		self.line_map['b_inte'].set_ydata(j * self.b_inte(u_arr))

		self.line_map['slope'].set_xdata([0.0, 50.0])
		self.line_map['slope'].set_ydata([0.0, 50.0 * self.slope])

		self.line_map['trace'].set_xdata(m['lon'])
		self.line_map['trace'].set_ydata(m['lat'])

		m['psi1'] = np.hstack(([math.nan,], (m['psi'][1:] - m['psi'][:-1]) / self.dg.dt))
		m['psi2'] = np.hstack(([math.nan,], (m['psi1'][1:] - m['psi1'][:-1]) / self.dg.dt))
		m['psi3'] = np.hstack(([math.nan,], (m['psi2'][1:] - m['psi2'][:-1]) / self.dg.dt))

		t_arr = np.array(m['t'])
		self.line_map['psidot_max'].set_xdata(t_arr)
		self.line_map['psidot_max'].set_ydata(np.ones_like(t_arr) * self.ct.psidot_limit(self.speed))

		self.line_map['psiddot_max'].set_xdata(t_arr)
		self.line_map['psiddot_max'].set_ydata(np.ones_like(t_arr) * self.ct.psiddot_limit(self.speed))

		for c in ['psi', 'psi1', 'psi2', 'psi3'] :
			self.line_map[c].set_xdata(t_arr)
			self.line_map[c].set_ydata(m[c])

		for k, v in self.axe_map.items() :
			v.relim()
			v.autoscale_view()

		self.fig.canvas.draw_idle()

	# def adapt(self, t) :
	# 	return np.clip(t**self.gamma, 0.0, 1.0)

	def b_func(self, t) :
		# u = self.adapt(t)
		return sum([self.c_lst[i]*t**i for i in range(self.n)])

	def b_inte(self, t) :
		# u = self.adapt(t)
		return sum([self.c_lst[i]*t**(i+1)/(i+1) for i in range(self.n)])

u = BezierPlot(10)
plt.show()
