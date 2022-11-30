#!/usr/bin/env python3

import collections
import math
from matplotlib.image import pil_to_array

import sympy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as pwd

from cc_pathlib import Path

from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider
from goto.aircraft.bspline import BSpline

""" l'idée est que de toute façon, si on defini psi comme une fonction de la distance à l'objectif alors, on défini aussi une variation de de cap, tout au long de la trajectoire.
reste à trouver une courbe qui ne fasse pas dépasser les limites de dérivées premières et seconde (et qui en plus est plutôt linéaire vers la fin)
"""


class ManualCurve() :

	period = 0.001

	def __init__(self, v_kt, psidot_max, phi_max, phidot_max) :
		self.line_map = dict()
		self.axe_map = dict()

		self.p_lst = [0.5,] * 24
		for i in range(16, 24) :
			self.p_lst[i] = 0.0
		print(self.p_lst)

		self.v_kt = v_kt
		self.v_ms = 1852 * v_kt / 3600

		self.psidot_max = psidot_max
		self.phi_max = phi_max
		self.phidot_max = phidot_max

		self.ct = CoordinatedTurn(psidot_max, phi_max)

		self.psidot_limit = self.ct.psidot_limit(self.v_ms)
		self.psiddot_limit = self.ct.psiddot_limit(self.v_ms)
		self.circle_radius = self.ct.circle_radius(self.v_ms)

		self.plot_setup()
		self.plot_update()

		plt.show()

	def run(self) :

		q_arr = np.cumsum([0.0,] + self.p_lst)
		y_arr = 90.0 * q_arr / q_arr[-1]
		
		self.dg = DummyGlider(lon=self.circle_radius * 1.667, vx=self.v_ms)

		x_arr = np.arange(25) / 16.0

		print("x_arr", len(x_arr), x_arr)
		print("y_arr", len(y_arr), y_arr)

		# plt.figure()
		# plt.plot(x_arr, y_arr)
		# plt.grid()

		b_arr = np.array([[x, y] for x, y in zip(x_arr, y_arr)])
		u = BSpline(b_arr, step=1.0)
		# print(b_arr, b_arr.shape)
		# print(u.timeline)

		plt.figure()
		plt.plot(b_arr[:,0], b_arr[:,1])
		c, d = u.trace()
		plt.plot(c[:,0], c[:,1])
		t_lst = list(range(25))
		plt.plot([1.5 * t / 24.0 for t in t_lst], [u.trace_curve(t)[1] for t in t_lst], '.')
		plt.grid()
		plt.show()

		for i in range(80000) :
			x = max(0.0, min(16 * self.dg.lon / self.circle_radius, 24.0))
			#psi = - np.interp(x, x_arr, y_arr)
			psi = - u.trace_curve(x)[1]
			print("psi", psi, x)
			self.dg.step_vx_psi(psi=psi)
			if self.dg.lon <= 0.25 and abs(self.dg.psi) <= 1.0 :
				break

		self.m_map = self.dg.freeze()

		return y_arr


	def plot_setup(self) :
		self.fig = plt.figure()

		self.axe_map["control"] = plt.subplot(2, 3, 1)
		plt.title("control")
		plt.grid()
		self.line_map['control_point'], = plt.plot([0.0,], 'o')
		self.line_map['control_line'], = plt.plot([0.0,], '.')

		self.axe_map["trace"] = plt.subplot(2, 3, 2)
		plt.title("trace")
		plt.grid()
		self.line_map['trace'], = plt.plot([0.0,])

		self.axe_map["psi"] = plt.subplot(3, 3, 3)
		plt.ylabel("psi")
		plt.grid()
		self.line_map['psi'], = plt.plot([0.0,])

		self.axe_map["psi1"] = plt.subplot(3, 3, 6)
		plt.ylabel("psi1")
		plt.grid()
		self.line_map['psi1'], = plt.plot([0.0,])
		self.line_map['psidot_max'], = plt.plot([0.0,], color='tab:red')

		self.axe_map["psi2"] = plt.subplot(3, 3, 9)
		plt.ylabel("psi2")
		plt.grid()
		self.line_map['psi2'], = plt.plot([0.0,])
		self.line_map['psiddot_max'], = plt.plot([0.0,], color='tab:red')

		self.a_lst = [
			self.fig.add_axes([0.05 + 0.045*(i % 12), 0.1 + ((24 - i - 1) // 12) * 0.2, 0.02, 0.12]) for i in range(24)
		]


		print(self.p_lst)
		self.s_lst = [
			pwd.Slider(ax=self.a_lst[i], label=f"P{i}", valmin=0.0, valmax=1.0, valinit=self.p_lst[i], orientation="vertical") for i in range(24)
		]

		for s in self.s_lst :
			s.on_changed(self.plot_update)


	def plot_update(self, value=0.0) :
		v_arr = np.array([s.val for s in self.s_lst])
		q_arr = np.cumsum(v_arr)
		p_arr = 90.0 * v_arr / q_arr[-1]

		self.p_lst = list(12.0 * p_arr / 90.0)

		for s, p in zip(self.s_lst, self.p_lst) :
			s.eventson = False
			s.set_val(p)
			s.eventson = True

		n_arr = np.arange(25)
		p_arr = self.run()

		m = self.m_map
		t_arr = np.array(m['t'])

		self.line_map['control_point'].set_xdata(n_arr)
		self.line_map['control_point'].set_ydata(p_arr)

		self.line_map['control_line'].set_xdata(16.0 * m['lon'] / self.circle_radius)
		self.line_map['control_line'].set_ydata(np.absolute(m['psi']))

		self.line_map['trace'].set_xdata(m['lon'])
		self.line_map['trace'].set_ydata(m['lat'])

		m['psi1'] = np.hstack(([math.nan,], (m['psi'][1:] - m['psi'][:-1]) / self.dg.dt))
		m['psi2'] = np.hstack(([math.nan,], (m['psi1'][1:] - m['psi1'][:-1]) / self.dg.dt))

		self.line_map['psidot_max'].set_xdata(t_arr)
		self.line_map['psidot_max'].set_ydata(np.ones_like(t_arr) * self.psidot_limit)

		self.line_map['psiddot_max'].set_xdata(t_arr)
		self.line_map['psiddot_max'].set_ydata(np.ones_like(t_arr) * self.psiddot_limit)

		for c in ['psi', 'psi1', 'psi2'] :
			self.line_map[c].set_xdata(t_arr)
			self.line_map[c].set_ydata(m[c])

		for k, v in self.axe_map.items() :
			v.relim()
			v.autoscale_view()

		self.fig.canvas.draw_idle()

	def test_direct(self, func, attr, d_ini=1000.0) :
		""" unlimited perfect response """
		x_lst, y_lst, psi_lst = list(), list(), list()
		x, y = 0, d_ini
		for i in range(50000) :
			x_lst.append(x)
			y_lst.append(y)
			psi = math.copysign( func(abs(y), ** attr), y )
			psi_lst.append(math.degrees(psi))
			x += self.v_ms * math.cos(-psi) * self.period
			y += self.v_ms * math.sin(-psi) * self.period

		self._save = {
			'x' : np.array(x_lst),
			'y' : np.array(y_lst),
			'psi' : np.array(psi_lst),
		}
		self._save['psid'] = (self['psi'][1:] - self['psi'][:-1]) / self.period
		self._save['psidd'] = (self['psid'][1:] - self['psid'][:-1]) / self.period

if __name__ == '__main__' :
	u = ManualCurve(100, 8.0, 24.0, 8.0)
