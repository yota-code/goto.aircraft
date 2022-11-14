#!/usr/bin/env python3


from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider

class PsiGo() :

	def __init__(self, speed=30.0) :
		self.speed = speed

		self.ct = CoordinatedTurn(12.0, 30.0)
		self.radius = self.ct.circle_radius(self.speed)

	def traj_linear(self, slope) :
		def traj(dist) :
			return min(slope / self.radius, 90.0)


	def run(self, func) :

		dist = 0.01
		self.dg = DummyGlider(lon=dist, vx=self.speed)

		while dist < 2 * self.radius :
			psi = func(dist)
			self.dg.step_vx_psi(psi=psi)
