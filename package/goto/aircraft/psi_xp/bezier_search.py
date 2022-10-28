#!/usr/bin/env python3


# search for an optimal solution

import math
import sys

import sympy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as pwd

from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider


class Bezier10() :
	def __init__(self, speed=30.0) :

		self.speed = speed
		self.slope = 1.0

		self.ct = CoordinatedTurn(12.0, 30.0)
		self.ct.display_info(self.speed)

		self.line_map = dict()
		self.axe_map = dict()
		self.score_map = dict()

	def detect(self) :
		p_arr = self.r_map['psi3']

		# first pass
		n = len(p_arr) // 32
		q_arr = np.convolve(p_arr, np.ones(n), mode='same') # strong smoothing

		a_arr = q_arr[:-2]
		u_arr = q_arr[1:-1]
		b_arr = q_arr[2:]

		x_arr = (u_arr - a_arr) * (b_arr - u_arr) # min/max detection
		m_arr = np.where(x_arr <= 0.0, np.ones_like(u_arr), np.zeros_like(u_arr))
		s_arr = np.where((b_arr - u_arr) <= 0.0, np.ones_like(u_arr), -np.ones_like(u_arr))

		r_arr = m_arr * s_arr
		self.r_map['peak_first_pass'] = r_arr

		r_lst = sorted([(i[0], 1) for i in np.argwhere(r_arr > 0.5)[-2:]] + [(i[0],-1) for i in np.argwhere(r_arr < -0.5)[-2:]])
		print(r_lst)

		p = len(p_arr)//512
		q_arr = np.convolve(p_arr, np.ones(p), mode='same')

		z_lst = list()
		for r, u in r_lst :
			z_arr = np.zeros_like(q_arr)
			e_arr = q_arr[r-n:r+n] * u
			self.r_map['peak_{r}_' + ('U' if u > 0 else 'D')] = r_arr
			z_arr[r-n:r+n] = e_arr - np.min(e_arr)
			z_lst.append((np.argmax(z_arr), u))
		print(z_lst)

		return z_lst

	def score(self) :
		if [i[1] for i in self.z_lst] == [-1, 1, -1, 1] :
			return (
				self.r_map['psi3'][self.z_lst[1][0]],
				self.r_map['psi2'][self.z_lst[1][0]]
			)
		else :
			return math.inf, math.inf


	def update(self, p_lst) :

		self.p_lst = [self.slope * i for i in p_lst]

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

	def evaluate(self, p_lst) :

		self.update(p_lst)

		self.r_map = self.run()
		self.z_lst = self.detect()

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


if __name__ == '__main__' :
	u = Bezier10()
	u.evaluate(1, 0.909, 0.0367, 0.416, 0.816, 0, 0.18, 0, 0, 0)