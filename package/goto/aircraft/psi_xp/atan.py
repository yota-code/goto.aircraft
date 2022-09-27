#!/usr/bin/env python3

import math

import numpy as np
import matplotlib.pyplot as plt

from goto.aircraft.coordinatedturn import CoordinatedTurn
from goto.aircraft.glider import DummyGlider

ct = CoordinatedTurn(12.0, 30.0)

spd = 1.0 # m.s-1
u = DummyGlider(lon=20.0, vx=spd)

def psi_atan(dst, spd, k1=32.0, k2=1.0) :
	return -k1 * math.atan(k2 * dst / spd)

for i in range(80000) :
	psi = psi_atan(u.lon, u.vx)
	u.step_vx_psi(psi=psi)
	if u.lon <= 0.25 and abs(u.psi) <= 1.0 :
		break

m = u.freeze()

print(m['t'])

plt.figure()
plt.subplot(2,2,1)
plt.plot(m['lon'], m['lat'])
plt.axis('equal')
plt.grid()
plt.subplot(2,2,2)
m['psidot'] = np.hstack(([math.nan,], (m['psi'][1:] - m['psi'][:-1]) / u.dt))
plt.plot(m['t'], m['psidot'])
plt.axhline(ct.psidot_deg(3600 * spd / 1852))
plt.ylabel('psidot')
plt.grid()
plt.subplot(2,2,3)
plt.plot(m['lon'], [psi_atan(dst, spd) for dst in m['lon']])
plt.xlabel('distance')
plt.ylabel('psi')
plt.grid()
plt.subplot(2,2,4)
m['psiddot'] = np.hstack(([math.nan,], (m['psidot'][1:] - m['psidot'][:-1]) / u.dt))
plt.plot(m['t'], m['psiddot'])
plt.ylabel('psiddot')
plt.grid()
plt.show()
