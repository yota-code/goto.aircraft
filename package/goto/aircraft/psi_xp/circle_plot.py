#!/usr/bin/env python3

import math

import sympy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as pwd

from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider


class CirclePlot() :

	def __init__(self) :

		self.ct = CoordinatedTurn(12.0, 30.0)
		self.speed = self.ct.transition_speed()

		self.line_map = dict()
		self.axe_map = dict()

		self.plot_setup()
		self.plot_update(0.0)

	def run(self, speed) :

		self.speed = speed

		circle_radius = self.ct.circle_radius(self.speed)

		self.dg = DummyGlider(lon=circle_radius * 1.1, vx=self.speed)

		for i in range(80000) :
			d = max(0.0, min(self.dg.lon / circle_radius, 1.0))
			psi = math.degrees(- math.acos(1.0 - d))
			if i == 0 :
				print(f"speed={self.speed} circle_radius={circle_radius}")
			self.dg.step_vx_psi(psi=psi)
			if self.dg.lon <= 0.25 and abs(self.dg.psi) <= 1.0 :
				break

		return self.dg.freeze()

	def plot_setup(self) :
		self.fig = plt.figure()

		self.axe_map["1.2.1"] = plt.subplot(1, 2, 1)
		self.line_map['trace'], = plt.plot([0.0,])
		plt.title("lat / lon")
		plt.grid()
		plt.axis('equal')

		self.axe_map["3.4.3"] = plt.subplot(4, 2, 2)
		self.line_map['psi'], = plt.plot([0.0,])
		plt.ylabel('psi')
		plt.grid()
		self.axe_map["3.4.7"] = plt.subplot(4, 2, 4)
		self.line_map['psi1'], = plt.plot([0.0,])
		plt.axhline(self.ct.psidot_limit(self.speed), color='tab:red')
		plt.ylabel('psidot (1)')
		plt.grid()
		self.axe_map["3.4.9"] = plt.subplot(4, 2, 6)
		self.line_map['psi2'], = plt.plot([0.0,])
		plt.axhline(self.ct.psiddot_limit(self.speed), color='tab:red')
		plt.ylabel('psiddot (2)')
		plt.grid()
		self.axe_map["3.4.12"] = plt.subplot(4, 2, 8)
		self.line_map['psi3'], = plt.plot([0.0,])
		plt.ylabel('psidddot (3)')
		plt.grid()

		self.a_lst = [
			self.fig.add_axes([0.1, 0.02, 0.8, 0.02]),
		]
		self.s_lst = [
			pwd.Slider(ax=a, label=f"speed", valmin=1.0, valmax=100.0, valinit=self.speed, orientation="horizontal") for i, a in enumerate(self.a_lst)
		]
		for s in self.s_lst :
			s.on_changed(self.plot_update)

	def plot_update(self, val) :

		m = self.run(* [s.val for s in self.s_lst])

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

		self.fig.canvas.draw_idle()


if __name__ == '__main__' :
	u = CirclePlot()
	plt.show()
