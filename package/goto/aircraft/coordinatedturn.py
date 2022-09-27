#!/usr/bin/env python3

import math

import numpy as np
import matplotlib.pyplot as plt

from cc_pathlib import Path

class CoordinatedTurn() :

	# prefered units : degrees and m/s

	g = 9.807
	e = 6371008.7714

	def __init__(self, psidot_max, phi_max, phidot_max=6.0) :
		
		self.psidot_max = psidot_max # in degrees
		self.phi_max = phi_max # in degrees
		self.phidot_max = phidot_max # in degrees

		self.vtr = self.transition_speed()

		print(f">>> psidot={self.psidot_max:0.3g}°/sec, phi={self.phi_max:0.3g}° transition={self.vtr:0.3g} m/s")

	def psidot_limit(self, spd) :
		return min( np.degrees(self.g * np.tan( np.radians(self.phi_max) ) / spd), self.psidot_max )

	def psiddot_limit(self, spd) :
		return self.psidot_limit(spd) * self.phidot_max / self.phi_max

	def transition_speed(self) :
		# return the transition speed in kt
		phi_rad = np.radians(self.phi_max)
		psidot_rad = np.radians(self.psidot_max)
		return self.g * np.tan(phi_rad) / psidot_rad

	def circle_radius(self, spd) :
		psidot_limit = self.psidot_limit(spd)
		psidot = np.where(spd < self.vtr, self.psidot_max, psidot_limit)
		t = math.tau / np.radians(psidot)
		return spd * t / math.tau

	def radius_to_curvature(self, radius) :
		aperture = radius / self.e
		curvature = 1 / ( np.tan(aperture) * self.e )
		return curvature

	def curvature_to_radius(self, curvature) :
		aperture = np.arctan(1 / (curvature * self.e))
		radius = aperture / self.e
		return radius

	def radius_to_max_speed(self, radius) :
		curvature = abs( self.radius_to_curvature(radius) )
		phi_rad = np.radians(self.phi_max)
		psidot_rad = np.radians(self.psidot_max)
		v_phi = math.sqrt( self.g * math.tan(phi_rad) / curvature )
		v_psidot = psidot_rad / curvature
		v_max = min(v_phi, v_psidot)
		print(f"radius = {radius:.6g}m -> curvature = {curvature:.4g} -> max speed = {v_max:.4g}m/s ({3600 * min(v_phi, v_psidot) / 1852:.6g} kt)")
		return v_max

	def max_speed_to_radius(self, v_kt) :
		psidot_deg = self.psidot_deg(v_kt)
		psidot_rad = np.radians(psidot_deg)
		v_ms = 1852 * v_kt / 3600
		curvature = psidot_rad / v_ms
		radius = self.radius_to_curvature(curvature)
		print(f"v_kt = {v_kt:.6g} kt ({v_ms:.4g}m/s) -> psidot = {psidot_deg:.4g}° -> curvature = {curvature:.4g}rad/m -> radius = {radius:.6g}m")
		return radius
			
def plot__turn_limit__() :

	ct_max = CoordinatedTurn(8.0, 20.0)
	ct_ffw = CoordinatedTurn(20.0 / 3.0, 50.0 / 3.0)

	sm_lst = np.linspace(0.0, ct_max.vtr)
	vm_lst = np.linspace(ct_max.vtr, 130)
	pm_lst = ct_max.psidot_deg(vm_lst)

	sf_lst = np.linspace(0.0, ct_ffw.vtr)
	vf_lst = np.linspace(ct_ffw.vtr, 130)
	pf_lst = ct_ffw.psidot_deg(vf_lst)

	plt.figure(figsize=(12, 8))
	plt.title(f"transition speed: max={ct_max.vtr:0.3f} / ffw={ct_ffw.vtr:0.3g} kt ")
	plt.plot(np.hstack((sm_lst, vm_lst)), np.hstack((np.ones_like(sm_lst) * ct_max.psidot, pm_lst)), color='tab:red', label="max") 
	plt.plot(np.hstack((sf_lst, vf_lst)), np.hstack((np.ones_like(sf_lst) * ct_ffw.psidot, pf_lst)), color='tab:green', label="ffw")
	plt.plot(np.hstack((sf_lst, vf_lst)), np.ones_like(np.hstack((sf_lst, vf_lst))) * 3.0, 'b--', label="3°/sec")
	plt.ylabel("psidot [°/dec]")
	plt.text(1.0, ct_max.psidot - 0.2, f"psidot_max={ct_max.psidot:0.3g}°/sec")
	plt.text(1.0, ct_ffw.psidot - 0.2, f"psidot_ffw={ct_ffw.psidot:0.3g}°/sec")
	v = (2*ct_max.vtr + 130) / 3
	plt.text(v - 4.0, ct_max.psidot_deg(v), f"phi_max={ct_max.phi:0.4g}°")
	v = (2*ct_ffw.vtr + 130) / 3
	plt.text(v - 4.0, ct_ffw.psidot_deg(v), f"phi_ffw={ct_ffw.phi:0.4g}°")

	plt.xlabel("air speed [kt]")
	plt.legend()
	plt.grid()
	plt.savefig("turn_limit.png")
	plt.show()

if __name__ == '__main__' :
	import sys
	import argparse

	parser = argparse.ArgumentParser(description='compute the min radius and the max speed admissibles')

	parser.add_argument('--radius_m', metavar='RADIUS', type=float, help='radius source')
	parser.add_argument('--speed_kt', metavar='SPEED', type=float, help='speed source')

	p = parser.parse_args()

	u = CoordinatedTurn(20.0 / 3.0, 50.0 / 3.0)
	if p.radius_m is not None :
		s = u.radius_to_max_speed(p.radius_m)
		r_check = u.max_speed_to_radius(3600 * s / 1852)
		print(r_check - p.radius_m)
	if p.speed_kt is not None :
		r = u.max_speed_to_radius(p.speed_kt)
		s_check = u.radius_to_max_speed(r)
		print(s_check - p.speed_kt)