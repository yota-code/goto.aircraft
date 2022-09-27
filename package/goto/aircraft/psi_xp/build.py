#!/usr/bin/env python3

"""
draw all possible trajectory for a given set of parameters
"""

import collections
import math
from ssl import ALERT_DESCRIPTION_HANDSHAKE_FAILURE
import sys

import numpy as np

from coordinated_turn import CoordinatedTurn

spd_kt = 100.0

phi_lim = 20.0
phidot_lim = 10.0
psidot_lim = 8.0
u = CoordinatedTurn(psidot_lim, phi_lim)

psidot_max = u.psidot_deg(spd_kt)
psiddot_max = phidot_lim / phi_lim * psidot_max

print(f"psidot_max = {psidot_max}")
print(f"psiddot_max = {psiddot_max}")

hsv_tpl = collections.namedtuple("HSV", ['hue', 'saturation', 'value'])
def hsv_cylindrical_interpolation(a:hsv_tpl, b:hsv_tpl, z:float) :
	# hue in radians
	# z is in 0.0, 1.0

	value = a.value * (1.0 - z) + b.value * z
	ax, ay = a.saturation * math.cos(a.hue), a.saturation * math.sin(a.hue)
	bx, by = b.saturation * math.cos(b.hue), b.saturation * math.sin(b.hue)

	x = ax * (1.0 - z) + bx * z
	y = ay * (1.0 - z) + by * z

	hue = math.atan2(y, x)
	saturation = math.sqrt(x*x + y*y)

	return hsv_tpl(hue, saturation, value)

def hsv_planar_interpolation(a:hsv_tpl, b:hsv_tpl, z:float) :
	# hue in radians
	# z is in 0.0, 1.0

	value = a.value * (1.0 - z) + b.value * z
	saturation = a.saturation * (1.0 - z) + b.saturation * z
	hue = ((a.hue * (1.0 - z) + b.hue * z) + math.tau) % math.tau

	return hsv_tpl(hue, saturation, value)





def integrate(cmd, jm, am, period) :
	a_lst = [0.0,]
	s_lst = [0.0,]
	t_lst = [0.0,]
	for i, j in enumerate(cmd) :
		a_lst.append(a_lst[-1] + j * jm * period )
		s_lst.append(s_lst[-1] + a_lst[-2] * period + j * jm * period * period * 0.5)
		t_lst.append((i +1) * period)
	return np.array(t_lst), np.array(a_lst), np.array(s_lst)

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
		x_lst.append(x_lst[-1] + math.cos(s_lst[-1]) * spd * period)
		y_lst.append(y_lst[-1] + math.sin(s_lst[-1]) * spd * period)

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

	return np.array(t_lst), np.array(a_lst), np.array(s_lst), np.array(x_lst), np.array(y_lst)

def cmd_iter(jm, am, period) :
	j_lst = list()
	while True :
		j = min( (am - (sum(j_lst) * jm * period)) / (jm * period), 1.0 )
		j_lst.append(j)
		j_arr = np.array(j_lst)
		yield np.hstack( (j_arr, - j_arr[::-1]) )

def plot_cmd(cmd, spd_kt, jm, am, period) :

	t_arr, a_arr, s_arr, x_arr, y_arr = navigate(j_arr, 1852 * spd_kt / 3600, jm, am, period)



if __name__ =='__main__' :
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import proj3d
	import matplotlib.colors as pcl

	jm = math.radians(psiddot_max)
	am = math.radians(psidot_max)
	period = 0.01

	print(f"am={am} jm={jm}")

	res_map = collections.defaultdict(set)

	m_arr = None
	for i, j_arr in enumerate( cmd_iter(jm, am, period) ) :

		if m_arr is None and j_arr[len(j_arr) // 2] == 0.0 :
			m_arr = p_arr
		
		t_arr, a_arr, s_arr, x_arr, y_arr = navigate(j_arr, 1852 * spd_kt / 3600, jm, am, period)

		for j, a, s, y in zip(j_arr, a_arr, s_arr, y_arr) :
			if s <= ( math.pi / 2.0 ) :
				res_map[round(j, 1)].add((y, s, a))
		
		if s_arr[-1] > 2.0 :
			break
		if i > 1500 :
			break

		p_arr = j_arr

	print(m_arr)

	left = hsv_tpl(0.0, 1.0, 1.0)
	right = hsv_tpl(math.pi, 1.0, 1.0)
	for z in range(11) :
		print(hsv_planar_interpolation(left, right, z / 10.0))

	fig = plt.figure()
	axe = fig.add_subplot(1, 1, 1, projection='3d')

	print(res_map.keys())

	for k, d in zip([-1.0, 0.0, 1.0], ['plum', 'royalblue', 'limegreen']) :
		#c = hsv_planar_interpolation(left, right, k / 2.0 + 0.5)
		#d = pcl.hsv_to_rgb((((c.hue / math.tau) + 1.0) % 1.0, c.saturation, c.value))
		p = np.array( sorted(res_map[k]) )

		axe.scatter(p[:,0], p[:,1], p[:,2], '.', color=d, alpha=0.1 if k == -1.0 else 1.0)

	t_arr, a_arr, s_arr, x_arr, y_arr = navigate(m_arr, 1852 * spd_kt / 3600, jm, am, period)

	# axe.scatter(y_arr, s_arr, a_arr, color='k', alpha=1.0)

	axe.set_xlabel("dist")
	axe.set_ylabel("angle")
	axe.set_zlabel("rate")
	plt.savefig("psi.png")
	plt.show()
	
	sys.exit(0)

	for j_lst in cmd_iter(am, jm, period) :
		u_lst = j_lst

