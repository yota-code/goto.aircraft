#!/usr/bin/env python3

import collections
import math

import sympy

import numpy as np
import matplotlib.pyplot as plt

from cc_pathlib import Path

from coordinated_turn import CoordinatedTurn

def func(x, y, c_lst, order) :
	o_lst = polynomial_basis(x, y, order)
	# print([c * o for c, o in zip(c_lst, o_lst)])
	return np.sum([c * o for c, o in zip(c_lst, o_lst)], axis=0)

def polynomial_basis(x, y, order, debug=False) :
	o_lst = list()

	try :
		ones = np.ones(len(x))
	except :
		ones = 1.0

	try :
		zeros = np.zeros(len(x))
	except :
		zeros = 0.0

	if debug :
		sx, sy = sympy.symbols('x y')
		ss = 0

	for i in range(order+1) :
		if i == 0 :
			o_lst.append(ones)
			if debug :
				ss += 1
		else :
			for j in range(i+1) :
				o_lst.append(
					(x**(i-j) if (0 < i-j) else ones) *
					(y**(j) if (0 < j) else ones)
				)

				if debug :
					ss += sx**(i-j) * sy**(j)
	if debug :
		print(ss)
	return np.array(o_lst)

class PhiFit() :
	period = 0.025
	def __init__(self, v_kt, psidot_max, phi_max, phidot_max) :
		self.v_kt = v_kt

		self.psidot_max = psidot_max
		self.phi_max = phi_max
		self.phidot_max = phidot_max

		self.ct = CoordinatedTurn(psidot_max, phi_max)

		self.am = self.ct.psidot_deg(v_kt)
		self.jm = self.am * phidot_max / phi_max 

		print(f"Am={self.am}\nJm={self.jm}")

	def cmd_maximal(self, tm) :
		# tm is the time at which the command is let at zero after the maximal rate has been reached
		j_lst = list()
		a_lst = [0.0,]
		s_lst = [0.0,]
		for i in range(n) :
			j = max(0.0, min(1.0, (self.am - a_lst[-1]) / (self.jm * self.period)))
			j_lst.append(j)
			a_lst.append(a_lst[-1] + j_lst[-1] * self.jm * self.period)
			s_lst.append(s_lst[-1] + a_lst[-2] * self.period + j_lst[-1] * self.jm * self.period * self.period * 0.5)

		cmd = j_lst + [(-j) for j in j_lst[::-1]]



	def generate_dataset(self) :

		# fig = plt.figure()
		# fig.suptitle(f"v_kt={self.v_kt}")
		# clr = plt.colormaps['plasma']

		m = collections.defaultdict(list)
		for b, cmd in enumerate( self.cmd__iter__() ) :
			if b == 150 :
				self.cmd_spe = cmd
			t_arr, j_arr, a_arr, s_arr, x_arr, y_arr = self.navigate(cmd)
			for n, (j, a, s, y) in enumerate(zip(j_arr, a_arr, s_arr, y_arr)) :
				if n == 0 :
					continue
				if s > 90.0 :
					continue
				m[j].append((y, s, a))

		w = {k : list() for k in m}
		for k in m :
			for y, s, a in m[k] :
				if y > 15.0 :
					w[k].append((y, s, a))

		for k in m :
			m[k] = np.array(sorted(w[k]))

		self.m = m

		Path(f"phifit_{self.v_kt}_{self.psidot_max}_{self.phi_max}.pickle.br").save(m)

		# for i, k in enumerate([-1.0, 0.0, 1.0]) :
		# 	axe = fig.add_subplot(2, 2, i+1, projection='3d')
		# 	axe.set_title(f"k={k}")
		# 	z = (k + 1.0) / 2.0
		# 	axe.scatter(m[k][:,0], m[k][:,1], m[k][:,2], '.', color=clr(z), alpha=self.period if k == -1.0 else 1.0)

		# axe = fig.add_subplot(2, 2, 4, projection='3d')
		# for k in m :
		# 	z = (k + 1.0) / 2.0
		# 	axe.scatter(m[k][:,0], m[k][:,1], m[k][:,2], '.', color=clr(z), alpha=self.period if k == -1.0 else 1.0)
		# axe.set_xlabel("dist")
		# axe.set_ylabel("angle")
		# axe.set_zlabel("rate")

		# # plt.savefig(f"phifit_{self.v_kt}_{self.psidot_max}_{self.phi_max}.png")
		# plt.close()

		# plt.figure()
		# axe = plt.axes(projection='3d')
		# z = (k + 1.0) / 2.0
		# k = -1.0
		# axe.scatter(m[k][:,1], m[k][:,2], m[k][:,0], color='b', alpha=0.01)
		# k = 1.0
		# axe.plot3D(m[k][:,1], m[k][:,2], m[k][:,0], '.', color='r')
		# plt.show()
		# plt.close()

		# self.fit_surface(m)

	def get_traj(self, cmd) :
		q = list()
		t_arr, j_arr, a_arr, s_arr, x_arr, y_arr = self.navigate(self.cmd_spe)
		for t, j, a, s, x, y in zip(t_arr, j_arr, a_arr, s_arr, x_arr, y_arr) :
			q.append((s, a, y))
		return q

	def fit_surface(self, order=5) :
		m = self.m

		k = -1.0
		x, y, z = m[k][:,1], m[k][:,2], m[k][:,0]

		# w = np.array([np.ones(len(z)), x, y, x**2, x * y, y**2])
		v = polynomial_basis(x, y, order, True)

		# plt.figure()
		# plt.subplot(3,2,1)
		# plt.plot(w[0,:] - v[0,:])
		# plt.subplot(3,2,2)
		# plt.plot(w[1,:] - v[1,:])
		# plt.plot(w[1,:] - x)
		# plt.subplot(3,2,3)
		# plt.plot(w[2,:] - v[2,:])
		# plt.subplot(3,2,4)
		# plt.plot(w[3,:] - v[3,:])
		# plt.subplot(3,2,5)
		# plt.plot(w[4,:] - v[4,:])
		# plt.subplot(3,2,6)
		# plt.plot(w[5,:] - v[5,:])

		# plt.show()

		# sys.exit()
		
		coef, resd, rank, sing = np.linalg.lstsq(v.T, z, rcond=None)

		p = list()
		for X in np.linspace(min(x), max(x), 32) :
			for Y in np.linspace(min(y), max(y), 32) :
				Z = func(X, Y, coef, order)
				p.append((X, Y, Z))
		p = np.array(p)

		q = np.array(self.get_traj(self.cmd_spe))
		
		fig = plt.figure()
		axe = fig.add_subplot(1, 2, 1, projection='3d', proj_type = 'ortho')
		axe.scatter(x, y, z, color='b', alpha=0.01)
		axe.scatter(p[:,0], p[:,1], p[:,2], color='r')
		axe.scatter(q[:,0], q[:,1], q[:,2], color='g')
		axe.set_xlabel("angle")
		axe.set_ylabel("rate")
		axe.set_zlabel("dist")
		axe = fig.add_subplot(1, 2, 2, projection='3d')
		axe.scatter(x, y, z - func(x, y, coef, order), alpha=0.01)
		#axe.scatter([0.0,], [0.0,], func(0.0, 0.0, coef, order), color='r')
		axe.set_xlabel("angle")
		axe.set_ylabel("rate")
		axe.set_zlabel("dist")
		plt.show()

	def cmd__iter__(self) :
		n = 1
		while True :
			j_lst = list()
			a_lst = [0.0,]
			s_lst = [0.0,]
			for i in range(n) :
				j = max(0.0, min(1.0, (self.am - a_lst[-1]) / (self.jm * self.period)))
				j_lst.append(j)
				a_lst.append(a_lst[-1] + j_lst[-1] * self.jm * self.period)
				s_lst.append(s_lst[-1] + a_lst[-2] * self.period + j_lst[-1] * self.jm * self.period * self.period * 0.5)

			cmd = j_lst + [(-j) for j in j_lst[::-1]]
			yield cmd
			if s_lst[-1] > 90.0 :
				break

			n += 1

	def navigate(self, cmd, debug=False) :

		v_ms = 1852 * self.v_kt / 3600

		j_lst = [0.0,]
		a_lst = [0.0,] # r
		s_lst = [0.0,] # psi
		t_lst = [0.0,] # time elapsed
		x_lst = [0.0,] # position x / y
		y_lst = [0.0,]
		for i, j in enumerate(cmd) :
			j_lst.append(j)
			t_lst.append((i+1) * self.period)
			a_lst.append(a_lst[-1] + j * self.jm * self.period )
			s_lst.append(s_lst[-1] + a_lst[-2] * self.period + j * self.jm * self.period * self.period * 0.5)
			x_lst.append(x_lst[-1] + math.cos(math.radians(s_lst[-1])) * v_ms * self.period)
			y_lst.append(y_lst[-1] + math.sin(math.radians(s_lst[-1])) * v_ms * self.period)

		t_arr, j_arr, a_arr, s_arr, x_arr, y_arr = np.array(t_lst), np.array(j_lst), np.array(a_lst), np.array(s_lst), np.array(x_lst), np.array(y_lst)

		if debug :
			plt.subplot(2,2,1)
			plt.title("J")
			plt.plot(cmd, '+--')
			plt.subplot(2,2,2)
			plt.title("A")
			plt.plot(t_arr, a_arr)
			plt.subplot(2,2,3)
			plt.title("S")
			plt.plot(t_arr, s_arr)
			plt.grid()
			plt.subplot(2,2,4)
			plt.title("traj")
			plt.plot(x_arr, y_arr)
			plt.axis('equal')
			plt.show()

		return t_arr, j_arr, a_arr, s_arr, x_arr, y_arr

if __name__ == "__main__" :


	# x = np.array([0.0, 1.0, 2.0])
	# y = np.array([0.0, 1.0, 2.0])
	# polynomial_basis(x, y, 2)

	# sys.exit()

	# for i in range(10, 100) :
	# 	u = PhiFit(i, 8.0, 20.0, 8.0)
	u = PhiFit(30.0, 8.0, 20.0, 8.0)
	u.generate_dataset()
	u.fit_surface()
	# for i, cmd in enumerate(u.cmd__iter__()) :
	# 	pass

	# u.navigate(cmd, True)