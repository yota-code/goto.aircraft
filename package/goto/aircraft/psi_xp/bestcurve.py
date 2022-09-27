#!/usr/bin/env python3

import collections
import math

import sympy

import numpy as np
import matplotlib.pyplot as plt

from cc_pathlib import Path

from aircraft import CoordinatedTurn

""" l'idée est que de toute façon, si on defini psi comme une fonction de la distance à l'objectif alors, on défini aussi une variation de de cap, tout au long de la trajectoire.
reste à trouver une courbe qui ne fasse pas dépasser les limites de dérivées premières et seconde (et qui en plus est plutôt linéaire vers la fin)
"""


class BestCurve() :
	period = 0.001
	def __init__(self, v_kt, psidot_max, phi_max, phidot_max) :
		self.v_kt = v_kt
		self.v_ms = 1852 * v_kt / 3600

		self.psidot_max = psidot_max
		self.phi_max = phi_max
		self.phidot_max = phidot_max

		self.ct = CoordinatedTurn(psidot_max, phi_max)

		self.am = self.ct.psidot_deg(v_kt)
		self.jm = self.am * phidot_max / phi_max

	def __getitem__(self, key) :
		return self._save[key]

	def pure_atan(self, d, k) :
		""" psi = atan(x)
		psid = 1 / (1+x**2)
		"""
		psi = math.atan(k * d / self.v_ms)
		return psi

	def weird_log(self, d, k) :
		psi = math.log((k * d / self.v_ms) + 1)
		return min(psi, math.pi / 2)

	def bezier_6(self, d, p0, p1, p2, p3, p4, p5) :
		"""
		import sympy
		p0, p1, p2, p3, p4, p5 = sympy.symbol('p0 p1 p2 p3 p4 p5')
		"""

		psi = (
			p0 * (1-d)**5 +
			p1 * 5 * (1-d)**4 * d + 
			p2 * 10 * (1-d)**3 * d**2 + 
			p3 * 10 * (1-d)**2 * d**3 + 
			p4 * 5 * (1-d)**1 * d**4 + 
			p5 * d**5
		) 
		return min(psi, math.pi / 2)

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
	u = BestCurve(100, 8.0, 20.0, 8.0)
	u.test_direct(u.pure_atan, {'k':0.3})
	v = BestCurve(100, 8.0, 20.0, 8.0)
	v.test_direct(u.weird_log, {'k':0.3})
	plt.subplot(2,2,1)
	plt.title("trajectory")
	plt.plot(u['x'], u['y'])
	plt.plot(v['x'], v['y'])
	plt.subplot(2,2,2)
	plt.title("psi [deg]")
	plt.plot(u['psi'])
	plt.plot(v['psi'])
	plt.subplot(2,2,3)
	plt.title("psid & psidd")
	plt.plot(u['psid'])
	plt.plot(u['psidd'])
	plt.plot(v['psid'])
	plt.plot(v['psidd'])
	plt.subplot(2,2,4)
	plt.title("psi = f(dist)")
	plt.plot(u['y'], u['psi'])
	plt.plot(v['y'], v['psi'])

	plt.show()	