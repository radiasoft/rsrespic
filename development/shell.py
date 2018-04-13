

import numpy as np
import matplotlib.pyplot as plt

import colormaps as cmaps

from solvers import * 
from particles import *
import analytic_models

sigma_x = 5.0
x_max = 20.
x = np.linspace(-x_max,x_max,100)
y = 0 * x

gaussian = analytic_models.uniform_gaussian(sigma_r = sigma_x, Q_0 = -2.0e-6)
r,phi = gaussian.compute_phi(x,y)

kv = analytic_models.kv(r_0 = sigma_x, Q_0 = -50.0e-6)
r,phi = kv.compute_phi(x,y)



## Solver parameters
lambda_0 = 60.0
n = 25

solver = field_solver_2D(lambda_x_0 = lambda_0, lambda_y_0 = lambda_0, 
	n_modes_x = n, n_modes_y = n )

## Particle distributions
particles = particles_2D(Q0 = 30000.0)
#particles.construct_gaussian(50000,sigma_x)
particles.construct_kv(50000, r_0 =sigma_x)


## Solver actions
solver.compute_phi(particles)
phi_grid, XX, YY = solver.compute_phi_mesh(xmax = x_max, ymax = x_max)


## post processing and  plotting 
phi_num = (phi_grid[:,50] - np.max(phi_grid[:,50])) #/ np.max(phi_grid[:,50] - phi_grid[0,50])
phi_a = phi  #/ np.max(phi - phi[0])

plt.figure()
plt.hexbin(particles.x,particles.y,cmap = cmaps.viridis) #plt.cm.RdBu)

plt.figure()
plt.plot(YY[:,50],phi_num, label = 'spectral solver')
plt.plot(x,phi_a, label = 'analytic model')
plt.legend(loc = 0)

plt.figure()
plt.plot(np.diff(phi_num) /np.diff(YY[:,50]))

plt.show()
