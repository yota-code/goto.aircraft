#!/usr/bin/env python3


# search for an optimal solution

import ast
import math
import sys

import sympy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as pwd

from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider

from cc_pathlib import Path


class Bezier10() :

	n = 10

	def __init__(self, speed=30.0) :

		self.speed = speed
		self.slope = 1.0

		self.ct = CoordinatedTurn(12.0, 30.0)
		self.ct.display_info(self.speed)

		self.line_map = dict()
		self.axe_map = dict()

		if Path("score.pson").is_file() :
			self.score_map = ast.literal_eval(Path("score.pson").read_text())
		else :
			self.score_map = dict()

	def detect(self) :

		m = self.r_map
		
		m['psi1'] = np.hstack(([math.nan,], (m['psi'][1:] - m['psi'][:-1]) / self.dg.dt))
		m['psi2'] = np.hstack(([math.nan,], (m['psi1'][1:] - m['psi1'][:-1]) / self.dg.dt))
		m['psi3'] = np.hstack(([math.nan,], (m['psi2'][1:] - m['psi2'][:-1]) / self.dg.dt))

		p_arr = m['psi3']

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
		m['peak_first_pass'] = r_arr

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
			return [
				self.r_map['psi3'][self.z_lst[1][0]],
				self.r_map['psi2'][self.z_lst[1][0]]
			]
		else :
			return [math.inf, math.inf]

	def plot(self, prefix) :
		m = self.r_map

		u_arr = np.linspace(0.0, 1.0, 512)
		b_arr = self.b_func(u_arr)

		j = 90.0 / self.final

		plt.figure(figsize=(18,12))

		plt.subplot(2, 3, 1)
		plt.plot(np.linspace(0.0, 1.0, self.n), self.p_lst, '.--')
		plt.plot(u_arr, b_arr)
		plt.grid()

		plt.subplot(2, 3, 4)
		plt.plot([0.0, 90.0, j * u_arr[-1]], [0.0, 90.0 * self.slope, 90.0])
		plt.plot(j * u_arr, j * self.b_inte(u_arr))
		plt.xlabel(j * u_arr[-1])
		plt.grid()

		plt.subplot(1, 3, 2)
		plt.plot(m['lon'], m['lat'])
		plt.title("lat / lon")
		plt.grid()
		plt.axis('equal')

		t_arr = np.array(m['t'])

		plt.subplot(4, 3, 3)
		plt.plot(t_arr, m['psi'])
		plt.ylabel("psi")
		plt.grid()

		plt.subplot(4, 3, 6)
		plt.plot(t_arr, m['psi1'])
		plt.axhline(self.ct.psidot_limit(self.speed), color='tab:red')
		plt.ylabel("psidot (1)")
		plt.grid()

		plt.subplot(4, 3, 9)
		plt.plot(t_arr, m['psi2'])
		plt.axhline(self.ct.psiddot_limit(self.speed), color='tab:red')
		plt.ylabel("psiddot (2)")
		plt.grid()

		plt.subplot(4, 3, 12)
		plt.plot([t_arr[i[0]] for i in self.z_lst], [m['psi3'][i[0]] for i in self.z_lst], 'o')
		plt.plot(t_arr, m['psi3'])
		plt.ylabel("psidddot (3)")
		plt.grid()

		plt.savefig(prefix + "_" + "_".join(f'{i:.03}' for i in self.p_lst) + '.png')
		plt.close()

	def update(self, p_lst) :

		self.p_lst = tuple(self.slope * i for i in p_lst)

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

	def evaluate(self, * p_lst) :

		self.update(p_lst)

		self.r_map = self.run()

		self.z_lst = self.detect()

		return self.score()

		# self.last_score = self.score()
		# self.score_map[self.p_lst] = self.last_score

		# Path("score.pson").write_text(repr(self.score_map))

		# print(p_lst, self.last_score)

		# self.plot()

	def run(self) :

		reference_radius = self.ct.circle_radius(30.0)
		circle_radius = self.ct.circle_radius(self.speed)

		j = 90.0 / self.final
		
		self.dg = DummyGlider(lon=circle_radius * j / reference_radius, vx=self.speed)

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

	def b_func(self, t) :
		# u = self.adapt(t)
		return sum([self.c_lst[i]*t**i for i in range(self.n)])

	def b_inte(self, t) :
		# u = self.adapt(t)
		return sum([self.c_lst[i]*t**(i+1)/(i+1) for i in range(self.n)])

if __name__ == '__main__' :

	def dst(score) :
		return math.sqrt(score[0]**2 + score[1]**2)

	def brb(p_lst) :
		return tuple([round(p, 3) for p in p_lst])

	u = Bezier10()

	if Path("score.pson").is_file() :
		score_map = ast.literal_eval(Path("score.pson").read_text())
	else :
		score_map = dict()

	b_val = (math.inf, math.inf)
	b_tup = None
	for p_tup in score_map :
		p_val = score_map[p_tup]
		if dst(p_val) < dst(b_val) :
			b_val = p_val
			b_tup = p_tup

	if b_tup is None :
		b_tup = brb([1.0, 1.0, 0.064, 0.445, 0.841, -0.011, 0.161, 0.0, 0.0, 0.0])
		
	b_val = u.evaluate(* b_tup)
	score_map[b_tup] = b_val
	u.plot(f"score({dst(b_val)})")

	sys.exit()

	for i in range(2**5) :

		d_lst = [0.0,] * 2 + [(1.0 if (i >> j) & 0x1 else -1.0) for j in range(5)] + [0.0,] * 3
		u_tup = brb([a + b*0.001 for a, b in zip(b_tup, d_lst)])
		if u_tup in score_map :
			continue

		u_val = u.evaluate(* u_tup)
		if dst(u_val) < dst(b_val) :
			score_map[u_tup] = u_val
			print(f"\x1b[32m{dst(u_val)} {u_val} {u_tup}\x1b[0m")
			u.plot(f"score({dst(u_val)})")
			
			b_val = u_val

	Path("score.pson").write_text(repr(score_map))

# 7.343454877134919e-05 [2.1316282072803006e-05, -7.027267656667391e-05] (1.0, 1.0, 0.06370000000000002, 0.445, 0.841, -0.011000000000000003, 0.16099999999999998, 0.0, 0.0, 0.0)
