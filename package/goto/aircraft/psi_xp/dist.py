#!/usr/bin/env python3

import math

import numpy as np
import matplotlib.pyplot as plt

from coordinated_turn import CoordinatedTurn

def navigate(cmd, spd, jm, am, period, debug=False) :
	a_lst = [0.0,] # r
	s_lst = [0.0,] # psi
	t_lst = [0.0,] # time elapsed
	x_lst = [0.0,] # position x / y
	y_lst = [0.0,]

	for i, j in enumerate(cmd) :
		a_lst.append(a_lst[-1] + j * jm * period )
		s_lst.append(s_lst[-1] + a_lst[-2] * period + j * jm * period * period * 0.5)
		t_lst.append((i +1) * period)
		x_lst.append(x_lst[-1] + math.cos(math.radians(s_lst[-1])) * spd * period)
		y_lst.append(y_lst[-1] + math.sin(math.radians(s_lst[-1])) * spd * period)

	t_arr, a_arr, s_arr, x_arr, y_arr = np.array(t_lst), np.array(a_lst), np.array(s_lst), np.array(x_lst), np.array(y_lst)

	if debug :
		plt.subplot(2,2,1)
		plt.title("J")
		plt.plot(cmd)
		plt.subplot(2,2,2)
		plt.title("A")
		plt.plot(t_arr, a_arr)
		plt.subplot(2,2,3)
		plt.title("S")
		plt.plot(t_arr, s_arr)
		plt.subplot(2,2,4)
		plt.title("traj")
		plt.plot(x_arr, y_arr)
		plt.axis('equal')
		plt.show()

	return t_arr, a_arr, s_arr, x_arr, y_arr

def distance_at_90degree(spd, am, jm) :
	t0 = am / jm
	s0 = am * t0 / 2.0
	s1 = 90.0 - ( 2.0 * s0 )
	t1 = s1 / am

	cmd = ( [jm,] * int(t0 * 1000) ) + ( [0.0,] * int(t1 * 1000) ) + ( [-jm,] * int(t0 * 1000) )

	print(t0, t1)
	print(s0, s1)

	t_arr, a_arr, s_arr, x_arr, y_arr = navigate(cmd, spd, jm, am, 0.001, debug=True)

	print([u[-1] for u in [t_arr, a_arr, s_arr, x_arr, y_arr]])

if __name__ == '__main__' :

	distance_at_90degree(5.0, 2.0, 1.0)