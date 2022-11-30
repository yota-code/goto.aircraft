#!/usr/bin/env python3

import sys
import numpy as np

def div2(n, d) :
	if d == 0.0 :
		return 0.0
	else :
		return n / d
		
def dist(v, norm=2) :
	return np.power(np.sum(np.power(v, norm)), 1/norm)
	
class BSpline(object) :
	def __init__(self, ctrl_pts, order=4, step=0.05):
		
		self.order = order
		self.step = step
				
		self.load_controls(ctrl_pts)
		
	def prepare_spline(self) :
		if self.k > self.n :
			self.k = self.n
		self.timeline = np.array(
			[0] * (self.k - 1) +
			list(range(self.n + 2 - self.k)) +
			[self.n + 1 - self.k] * (self.k-1), dtype=np.float64)
		self.t = []
		
	def load_controls(self, ctrl_pts) :
		self.calculated_curve = None
		self.n = 0
		self.k = self.order
		self.timeline = None
		
		self.dimension = len(ctrl_pts[0])

		self.curve = np.empty((0, self.dimension), dtype=np.float64)
		self.derivate = np.empty((0, self.dimension), dtype=np.float64)
		self.controls = np.empty((0, self.dimension), dtype=np.float64)
		
		for p in ctrl_pts :
			self.add_control(p)
			
		self.prepare_spline()
		
	def add_control(self, ctrl_pt) :
		self.controls = np.vstack((self.controls, ctrl_pt))
		self.n += 1
			
	def trace(self) :
		for u in range(len(self.timeline) - 1) :
			self.t.append(self.timeline[u])
			while self.t[-1] < self.timeline[u+1] :
				self.t.append(self.walk(self.t[-1]))
		self.t.append(self.walk(self.timeline[u+1]))
		return self.curve, self.derivate
		
	def walk(self, t) :
		_c = self.trace_curve(t)
		print(t, _c)
		self.curve = np.vstack((self.curve, _c))
		
		_d = self.trace_derivate(t)
		self.derivate = np.vstack((self.derivate, _d))
		
		return t + self.step / max(dist(_d), self.step)
		
	def trace_curve(self, t) :
		if self.timeline[-1] <= t :
			il = self.n
			t = self.timeline[-1]
		else :
			il = 0
			while not ( self.timeline[il] <= t < self.timeline[il+1] ) :
				il += 1
			il += 1
			
		self.N0 = np.zeros((self.n+1,self.k), dtype=np.float64)
		
		self.N0[il-1,0] = 1
		ti = self.timeline
		for k in range(2, self.k+1) :
			for i in range(il - k, il):
				u =  self.N0[i, k-2] * div2(( t - ti[i]),(ti[i+k-1] - ti[i]))
				u += self.N0[i+1, k-2] * div2((ti[i+k] - t),(ti[i+k] - ti[i+1]))
				self.N0[i,k-1] = u
				
		self.calculated_curve = t
		self.il = il
		
		return np.dot(self.N0[:-1,-1], self.controls)
				
	def trace_derivate(self, t) :
		if self.calculated_curve != t :
			self.curve(t)
		il = self.il

		self.N1 = np.zeros((self.n+1,self.k), dtype=np.float64)
		
		ti = self.timeline
		for k in range(2,self.k+1) :
			for i in range(il - k, il):
				u =  self.N1[i  , k-2] * div2((  t - ti[i]), (ti[i+k-1] - ti[i]  ))
				u += self.N0[i  , k-2] * div2((          1), (ti[i+k-1] - ti[i]  ))
				u += self.N1[i+1, k-2] * div2((ti[i+k] - t), (ti[i+k]   - ti[i+1]))
				u += self.N0[i+1, k-2] * div2((         -1), (ti[i+k]   - ti[i+1]))
				self.N1[i,k-1] = u

		return np.dot(self.N1[:-1,-1], self.controls)
		
if __name__ == '__main__' :

	u = np.array([
		[-1.0, -1.0],
		[0.0, -1.0],
		[0.0, 1.0],
		[1.0, 1.0],
		[2.0, -1.0],
		[2.0, 1.0],
	])

	u = np.array([
		[0.0,],
		[1.0,],
		[2.0,],
		[4.0,],
		[5.0,],
		[6.0,],
	])

	u = np.array([
		[0.0, 0.0,],
		[1.0, 0.0,],
		[2.0, 0.0,],
		[3.0, 1.0,],
		[4.0, 1.0,],
		[5.0, 1.0,],
	])


	
	v = BSpline(u)
	c, d = v.trace()
	#t_arr = np.linspace(2.0, 3.0)

	#y_lst = [v.trace_curve(t) for t in t_arr]

	import matplotlib.pyplot as plt

	plt.plot(u[:,0], u[:,1])
	plt.plot(c[:,0], c[:,1])
	#plt.plot(t_arr, y_lst)
	plt.grid()
	plt.show()
