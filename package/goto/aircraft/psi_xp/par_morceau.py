import numpy as np
import matplotlib.pyplot as plt

import sympy

x, b, e = sympy.symbols('x b e')

U0, V0, W0, U1, V1, W1 = sympy.symbols(' '.join([f'{j}{i}' for i in range(2) for j in 'UVW']))

p_lst = sympy.symbols(' '.join(f'p{i}' for i in range(6)))
x_lst = [x**i for i in range(6)]

u = sum((p * x) for p, x in zip(p_lst, x_lst))
v = u.diff('x')
w = v.diff('x')

x0 = {'x': b}
x1 = {'x': e}

r = sympy.solve([u.subs(x0) - U0, u.subs(x1) - U1, v.subs(x0) - V0, v.subs(x1) - V1, w.subs(x0) - W0, w.subs(x1) - W1], ['p0', 'p1', 'p2', 'p3', 'p4', 'p5'])

# print(r)

z = {'b': 1, 'e': 8, 'U0': 3, 'U1': 7, 'V0': 1, 'V1': -1, 'W0':1, 'W1': 0}

u_sol = u.subs(r).subs(z)

t_arr = np.linspace(z['b'], z['e'])
y_arr = np.array([u_sol.subs({'x': t}) for t in t_arr])

plt.plot(t_arr, y_arr)
plt.grid()
plt.show()