#!/usr/bin/env python3

import collections
import math

import numpy as np

from cc_pathlib import Path

""" the earth is flat
the gravity is equal to 9.807 everywhere
all data are expressed in SI
the model is only valid for perfectly coordinated turns
"""

earth_gravity = 9.807

class DummyGlider() :

	dt_ms = 1

	rec_lst = ['t', 'lat', 'lon', 'alt', 'psi', 'vx', 'vz']

	def __init__(self, lat=0.0, lon=0.0, alt=100.0, psi=0.0, vx=40.0, vz=0.0) :

		self.t = 0
		self.dt = self.dt_ms / 1000.0

		self.lat = lat # in meter
		self.lon = lon
		self.alt = alt
		self.psi = psi
		self.vx = vx
		self.vz = vz

		self.rec = collections.defaultdict(list)

	def record(self) :
		for k in self.rec_lst :
			self.rec[k].append( getattr(self, k) )
		self.rec['t'][-1] *= self.dt
		self.t += self.dt_ms

	def freeze(self) :
		return {
			k: np.array(v) for k, v in self.rec.items()
		}

	def save(self, pth) :

		pth = Path(pth).with_suffix('.tsv')

		stack = [
			self.rec_lst,
		]
		for i in range(len(self.rec['t'])) :
			stack.append(
				[ self.rec[k][i] for k in self.rec_lst ]
			)

		pth.save(stack)

	def step(self, vx, rol, vz=0.0) :

		if vx is None :
			vx = self.vx
		if vz is None :
			vz = self.vz

		self.t += self.dt_ms

		psi_dot = earth_gravity * math.tan(rol) / self.vx

		self.psi += psi_dot * self.dt

		self.vx = vx

		self.lat += self.vx * math.cos(self.psi) * self.dt
		self.lon += self.vx * math.sin(self.psi) * self.dt

		self.vz = vz

		self.alt += self.vz * self.dt

		self.record()

	def step_vx_psi(self, vx=None, psi=None, vz=None) :

		if vz is not None :
			self.vz = vz
		self.alt += self.vz * self.dt

		if vx is not None :
			self.vx = vx
		if psi is not None :
			self.psi = psi

		self.lat += self.vx * math.cos(math.radians(self.psi)) * self.dt
		self.lon += self.vx * math.sin(math.radians(self.psi)) * self.dt

		self.record()
